"""Part 1's working figures: fig0 to fig6, E_regression and F_maps (retired).

PNG at 200 dpi beside the tables of adult_v2/ish_analysis/beyond/, for checking a
step as it ran: fig0_structures to fig3_residual drawn by run_beyond_density.py
(the structures of the fit, the ceiling, what each model predicts, the leftover's
replication and where it is largest), fig4_controls to fig6_readings by
run_beyond_controls.py (controls A to D, E and F, G) and E_regression and F_maps by
run_beyond_regression.py (the fit with its diagnostic, the three maps on the
brain). The guided figures 03, 03s1, 03s2 and 04 show the same numbers for a
reader, and the steps print them. The run scripts drew these from what the steps'
main() returned, which no longer returns them.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from sepmap.adult.plotting import saved_working
from sepmap.plotting import DARK_BLUE, DARK_GREY, MID_GREY, RED, paint, tidy

# ===== Working figures of analysis 4 (adult.beyond_density, beyond_controls) =====

# the readings of control G as the working figure names them: the map's own zref, the
# same with the young-against-adult tables' 17-brain reference, and the stored
# readings
READING_LABELS = {
    "zref": "zref\n(declared reference)",
    "zref_stored": "zref\n(17-brain reference)",
    "cref": "cref",
    "subref": "subref",
    "ratio": "nano / auto",
    "sepratio": "nano / SEP",
}


def plot_structures_used(rows: pd.DataFrame, save: Path | None = None) -> plt.Figure:
    """fig0: the structures used per division, and those left out per reason.

    `rows` is structures_used.csv (beyond_density.structure_rows).
    """
    fig, axes = plt.subplots(
        1, 2, figsize=(11.0, 3.9), gridspec_kw=dict(width_ratios=[1, 1.25])
    )
    counts = rows[rows["used"]].groupby("division").size().sort_values(ascending=False)
    axes[0].bar(
        range(len(counts)), counts.to_numpy(), color="0.65", edgecolor="0.25", lw=0.5
    )
    axes[0].set_xticks(range(len(counts)))
    axes[0].set_xticklabels(counts.index, fontsize=7, rotation=45, ha="right")
    axes[0].set_ylabel("structures used", fontsize=8)
    axes[0].set_title(
        f"what the fit runs on\n{int(rows['used'].sum())} grey-matter structures",
        fontsize=9,
    )
    tidy(axes[0])

    # the reasons, the declared set's own reasons pooled into one
    why = rows.loc[~rows["used"], "reason"].map(
        lambda r: "not in the declared set" if r.startswith("not in the declared") else r
    )
    why_counts = why.value_counts().sort_values()
    axes[1].barh(
        range(len(why_counts)),
        why_counts.to_numpy(),
        color=RED,
        edgecolor="0.25",
        linewidth=0.5,
        alpha=0.85,
    )
    axes[1].set_yticks(range(len(why_counts)))
    axes[1].set_yticklabels([w[:52] for w in why_counts.index], fontsize=7)
    axes[1].set_xlabel("structures left out", fontsize=8)
    axes[1].set_title(
        "and what it leaves out\ndecided by rule, before any fitting", fontsize=9
    )
    tidy(axes[1])
    fig.tight_layout()
    return saved_working(fig, save)


def plot_ceiling(
    agreement: list[float],
    explainable: float,
    reproducible: bool,
    save: Path | None = None,
) -> plt.Figure:
    """fig1: the two half-cohort maps' agreement over every split, and the ceiling.

    `reproducible` says whether the mean agreement replicates
    (beyond_density.replicates), which the title follows.
    """
    half = float(np.mean(agreement))
    fig, ax = plt.subplots(figsize=(5.6, 3.9))
    ax.hist(agreement, bins=25, color="0.7", edgecolor="0.35", linewidth=0.4)
    ax.axvline(half, color=RED, lw=1.8)
    ax.set_xlabel("Spearman between the two half-cohort maps", fontsize=8)
    ax.set_ylabel(f"splits of ten adults ({len(agreement)})", fontsize=8)
    verdict = "is reproducible" if reproducible else "does not reproduce"
    ax.set_title(
        f"Step 1. the map {verdict}\n"
        f"half-cohorts agree at {half:.3f}; ceiling (Spearman-Brown) {explainable:.3f}",
        fontsize=9,
    )
    tidy(ax)
    fig.tight_layout()
    return saved_working(fig, save)


def plot_covariates(
    table: pd.DataFrame, explainable: float, save: Path | None = None
) -> plt.Figure:
    """fig2: each model's cross-validated R2 against the ceiling.

    `table` is variance_partition.csv (beyond_density.partition_table).
    """
    fig, ax = plt.subplots(figsize=(7.8, 4.5))
    ax.barh(np.arange(len(table)), table["cv_r2"], color="0.65", edgecolor="0.25", lw=0.5)
    ax.axvline(explainable, color=RED, lw=1.8)
    ax.annotate(
        f"ceiling {explainable:.3f}\n(the map's own reliability)",
        (explainable, len(table) - 0.4),
        color=RED,
        fontsize=7.5,
        ha="right",
        va="top",
        xytext=(-6, 0),
        textcoords="offset points",
    )
    ax.set_yticks(np.arange(len(table)))
    ax.set_yticklabels(table["model"], fontsize=7.5)
    ax.set_xlim(min(0.0, float(table["cv_r2"].min()) - 0.02), 1.02)
    ax.set_xlabel(
        "variance of the map predicted on held-out structures (cross-validated R2)",
        fontsize=8,
    )
    quoted = float(table.set_index("key").loc["model", "cv_r2"])
    ax.set_title(
        f"Step 2. what Gria1 and synapse density predict\n"
        f"the main model: {quoted / explainable:.1%} of the ceiling, "
        f"{1 - quoted / explainable:.1%} left",
        fontsize=9,
    )
    tidy(ax)
    fig.tight_layout()
    return saved_working(fig, save)


def plot_residual(
    map_agreement: list[float],
    agreement: list[float],
    implied: float,
    table: pd.DataFrame,
    save: Path | None = None,
) -> plt.Figure:
    """fig3: the half-cohort agreements, and the structures most off prediction.

    `map_agreement` and `agreement` are the maps' and the leftovers' agreement per
    split, `implied` the replication the ceiling and the fit alone imply, `table`
    residual_by_structure.csv.
    """
    fig, axes = plt.subplots(
        1, 2, figsize=(11.8, 4.5), gridspec_kw=dict(width_ratios=[1, 1.15])
    )
    bins = np.linspace(
        min(min(map_agreement), min(agreement)) - 0.005,
        max(max(map_agreement), max(agreement)) + 0.005,
        30,
    )
    axes[0].hist(
        map_agreement,
        bins=bins,
        color="0.6",
        edgecolor="0.3",
        lw=0.3,
        alpha=0.85,
        label="the map itself",
    )
    axes[0].hist(
        agreement,
        bins=bins,
        color=RED,
        edgecolor="0.3",
        lw=0.3,
        alpha=0.7,
        label="what is left of it",
    )
    axes[0].axvline(implied, color="0.1", lw=1, ls=(0, (3, 2)))
    axes[0].set_xlabel("half-cohort against half-cohort (Spearman)", fontsize=8)
    axes[0].set_ylabel(f"splits of ten adults ({len(agreement)})", fontsize=8)
    axes[0].legend(fontsize=7.5, frameon=False, loc="upper left")
    axes[0].set_title(
        f"Step 3. map {np.mean(map_agreement):.3f}, leftover {np.mean(agreement):.3f}\n"
        f"dashed: {implied:.3f}, what the ceiling and the fit alone imply",
        fontsize=9,
    )
    tidy(axes[0])

    show = pd.concat([table.head(8), table.tail(6)])
    axes[1].barh(
        np.arange(len(show)),
        show["residual"],
        color=[RED if v > 0 else DARK_BLUE for v in show["residual"]],
        edgecolor="0.25",
        linewidth=0.4,
    )
    axes[1].set_yticks(np.arange(len(show)))
    axes[1].set_yticklabels([s[:36] for s in show["structure"]], fontsize=6.8)
    axes[1].invert_yaxis()
    axes[1].axvline(0, color="0.3", lw=0.7)
    axes[1].set_xlabel("nano rank minus predicted rank", fontsize=8)
    axes[1].set_title("Step 4. where it is largest", fontsize=9)
    tidy(axes[1])
    fig.tight_layout()
    return saved_working(fig, save)


def plot_artefacts(
    leftover: np.ndarray,
    xyz: np.ndarray,
    sizes: np.ndarray,
    pairs: list[float],
    naive: np.ndarray,
    rws: np.ndarray,
    passed: dict[str, bool],
    save: Path | None = None,
) -> plt.Figure:
    """fig4: controls A to D, position, size, pairs of adults, naive against RWS.

    `xyz` holds each structure's centroid (AP, DV, ML) in mm, `passed` each
    control's verdict by letter; each title says what its verdict says.
    """
    fig, axes = plt.subplots(1, 4, figsize=(15.5, 3.9))

    # A: the leftover against anterior-posterior position
    axes[0].scatter(
        xyz[:, 0], leftover, s=12, facecolor="0.6", edgecolor="0.25", linewidth=0.3
    )
    axes[0].axhline(0, color="0.85", lw=0.7)
    axes[0].set_xlabel("structure centroid, anterior-posterior (mm)", fontsize=8)
    axes[0].set_ylabel("residual (ranks)", fontsize=8)
    if passed["A"]:
        title = "A. not a front-to-back gradient"
    else:
        title = "A. a gradient could explain it"
    axes[0].set_title(title, fontsize=9)
    tidy(axes[0])

    # B: the leftover against structure volume
    good = np.isfinite(sizes)
    axes[1].scatter(
        np.log10(sizes[good]),
        leftover[good],
        s=12,
        facecolor="0.6",
        edgecolor="0.25",
        linewidth=0.3,
    )
    axes[1].axhline(0, color="0.85", lw=0.7)
    axes[1].set_xlabel("log10 structure volume (20 um voxels)", fontsize=8)
    axes[1].set_ylabel("residual (ranks)", fontsize=8)
    title = "B. not small-structure noise" if passed["B"] else "B. size may drive it"
    axes[1].set_title(title, fontsize=9)
    tidy(axes[1])

    # C: the agreement of every pair of adults, the median in red
    axes[2].hist(pairs, bins=20, color="0.7", edgecolor="0.35", linewidth=0.4)
    axes[2].axvline(float(np.median(pairs)), color=RED, lw=1.8)
    axes[2].set_xlabel("leftover of one adult against another (Spearman)", fontsize=8)
    axes[2].set_ylabel("pairs of adults", fontsize=8)
    title = "C. every animal shows it" if passed["C"] else "C. one animal may carry it"
    axes[2].set_title(title, fontsize=9)
    tidy(axes[2])

    # D: the naive leftover against the RWS one, with the identity line
    axes[3].scatter(naive, rws, s=12, facecolor="0.6", edgecolor="0.25", linewidth=0.3)
    lim = [min(naive.min(), rws.min()) - 3, max(naive.max(), rws.max()) + 3]
    axes[3].plot(lim, lim, color="0.75", ls="--", lw=0.8)
    axes[3].set_xlabel("leftover, five naive adults", fontsize=8)
    axes[3].set_ylabel("leftover, five RWS adults", fontsize=8)
    title = (
        "D. not the whisker manipulation" if passed["D"] else "D. naive and RWS disagree"
    )
    axes[3].set_title(f"{title}\nrho {spearmanr(naive, rws).statistic:+.2f}", fontsize=9)
    tidy(axes[3])

    not_ruled_out = [k for k in "ABCD" if not passed[k]]
    if not_ruled_out:
        title = (
            "Four ways the leftover could be an artefact; not ruled out: "
            + ", ".join(not_ruled_out)
        )
    else:
        title = "Four ways the leftover could be an artefact, and is not"
    fig.suptitle(title, fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    return saved_working(fig, save)


def plot_model_space(
    curve: pd.DataFrame,
    k_most: int,
    explainable: float,
    scores: tuple[float, float, float],
    passed: dict[str, bool],
    save: Path | None = None,
) -> plt.Figure:
    """fig5: controls F and E, the gene-space curve and curving the model.

    `curve` is gene_space.csv, `k_most` the number of components picked most often,
    `scores` the held-out R2 of the model straight, curved to cubes and to fifth
    powers.
    """
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    ks = curve["n_components"]
    axes[0].plot(ks, curve["r2"], color="0.6", lw=1.5, label="fitted")
    axes[0].plot(ks, curve["cv_r2"], color=RED, lw=1.8, label="cross-validated")
    axes[0].axhline(explainable, color="0.3", ls="--", lw=1.2)
    axes[0].annotate(
        "ceiling", (ks.iloc[-1], explainable), fontsize=7.5, ha="right", va="bottom"
    )
    axes[0].axvline(k_most, color="0.4", ls=":", lw=1.0)
    axes[0].set_xlabel(
        "components of the expression of the genes measured in every structure",
        fontsize=8,
    )
    axes[0].set_ylabel("variance of the map explained", fontsize=8)
    axes[0].legend(fontsize=7.5, frameon=False)
    if passed["F"]:
        title = "F. even the whole gene table falls short"
    else:
        title = "F. the whole gene table accounts for the map"
    axes[0].set_title(
        f"{title}\nthe gap between the two lines is overfitting", fontsize=9
    )
    tidy(axes[0])

    axes[1].bar(
        [0, 1, 2], scores, color=[RED, "0.65", "0.8"], edgecolor="0.25", linewidth=0.5
    )
    axes[1].set_xticks([0, 1, 2])
    axes[1].set_xticklabels(
        ["the model\n(straight)", "curved\n(squares, cubes)", "to fifth powers"],
        fontsize=8,
    )
    axes[1].set_ylabel("cross-validated R2", fontsize=8)
    if passed["E"]:
        title = "E. and curving it buys nothing"
    else:
        title = "E. and curving it pays: read the check row curved"
    axes[1].set_title(title, fontsize=9)
    tidy(axes[1])
    fig.tight_layout()
    return saved_working(fig, save)


def plot_readings(
    table: pd.DataFrame, passed: bool, save: Path | None = None
) -> plt.Figure:
    """fig6: control G, per reading the map's and the leftover's replication, and R2.

    `table` is readings.csv.
    """
    fig, ax = plt.subplots(figsize=(8.0, 4.0))
    x = np.arange(len(table))
    ax.bar(
        x - 0.22,
        table["map_replication"],
        width=0.2,
        color="0.55",
        edgecolor="0.25",
        linewidth=0.4,
        label="the map replicates",
    )
    ax.bar(
        x,
        table["leftover_replication"],
        width=0.2,
        color=RED,
        edgecolor="0.25",
        linewidth=0.4,
        label="the leftover replicates",
    )
    ax.bar(
        x + 0.22,
        table["r2"],
        width=0.2,
        color="0.8",
        edgecolor="0.25",
        linewidth=0.4,
        label="the model explains (R2)",
    )
    ax.set_xticks(x)
    ax.set_xticklabels([READING_LABELS[r] for r in table["reading"]], fontsize=7.5)
    ax.set_ylim(0, 1.05)
    ax.legend(fontsize=7.5, frameon=False, loc="lower right")
    if passed:
        title = "G. the same picture under every reading, not just zref"
    else:
        title = "G. not the same picture under every reading"
    ax.set_title(title, fontsize=9)
    tidy(ax)
    fig.tight_layout()
    return saved_working(fig, save)


# ===== Working figures of the regression (adult.beyond_regression) =====


def draw_fit(
    ax: plt.Axes,
    observed: np.ndarray,
    predicted: np.ndarray,
    leftover: np.ndarray,
    structures: list[str],
    fitted_r2: float,
    cv: float,
) -> None:
    """Draw observed against predicted, the five largest residuals named."""
    lim = [
        min(observed.min(), predicted.min()) - 4,
        max(observed.max(), predicted.max()) + 4,
    ]
    ax.plot(lim, lim, color=MID_GREY, ls="--", lw=1.0, zorder=1)
    ax.scatter(
        predicted,
        observed,
        s=16,
        facecolor=DARK_GREY,
        edgecolor="0.2",
        linewidth=0.3,
        zorder=2,
    )
    worst = np.argsort(-np.abs(leftover))[:5]
    for i in worst:
        ax.annotate(
            structures[i][:24],
            (predicted[i], observed[i]),
            fontsize=6.5,
            color=RED,
            xytext=(4, 2),
            textcoords="offset points",
        )
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel("predicted from Gria1 and synapse density (rank)", fontsize=8.5)
    ax.set_ylabel("nano map (rank)", fontsize=8.5)
    ax.set_title(f"the fit\nR2 {fitted_r2:.2f} fitted, {cv:.2f} predicted", fontsize=9.5)
    tidy(ax)


def draw_diagnostic(
    ax: plt.Axes,
    predicted: np.ndarray,
    leftover: np.ndarray,
    misspecified: bool,
    p_max: float,
) -> None:
    """Draw the residual against the prediction, the standard diagnostic.

    The title follows `misspecified`, whether a tilt or a fan reached `p_max`.
    """
    ax.axhline(0, color=MID_GREY, lw=1.0)
    ax.scatter(
        predicted,
        leftover,
        s=16,
        facecolor=DARK_GREY,
        edgecolor="0.2",
        linewidth=0.3,
        zorder=2,
    )
    ax.set_xlabel("predicted (rank)", fontsize=8.5)
    ax.set_ylabel("residual (ranks)", fontsize=8.5)
    if misspecified:
        title = (
            f"the diagnostic\na tilt or a fan (p < {p_max:g}), so the model\n"
            "may be mis-specified"
        )
    else:
        title = (
            "the diagnostic\nno tilt and no fan, so the model is\n"
            "incomplete rather than mis-specified"
        )
    ax.set_title(title, fontsize=9.5)
    tidy(ax)


def draw_leftover(ax: plt.Axes, leftover: np.ndarray) -> None:
    """Draw the distribution of the residuals, the leftover."""
    ax.hist(leftover, bins=26, color=MID_GREY, edgecolor="0.3", linewidth=0.4)
    ax.axvline(0, color=DARK_GREY, lw=1.4)
    ax.set_xlabel("residual (ranks)", fontsize=8.5)
    ax.set_ylabel("structures", fontsize=8.5)
    ax.set_title(
        f"the leftover\nspread {leftover.std():.1f} ranks over {len(leftover)} "
        "structures",
        fontsize=9.5,
    )
    tidy(ax)


def plot_regression(
    observed: np.ndarray,
    predicted: np.ndarray,
    leftover: np.ndarray,
    structures: list[str],
    fitted_r2: float,
    cv: float,
    misspecified: bool,
    p_max: float,
    save: Path | None = None,
) -> plt.Figure:
    """E_regression: observed against predicted, the diagnostic, the residuals.

    `misspecified` says whether a tilt or a fan of the residuals reached `p_max`
    (beyond_regression.diagnostic).
    """
    fig, axes = plt.subplots(1, 3, figsize=(13.4, 4.3))
    draw_fit(axes[0], observed, predicted, leftover, structures, fitted_r2, cv)
    draw_diagnostic(axes[1], predicted, leftover, misspecified, p_max)
    draw_leftover(axes[2], leftover)
    fig.suptitle(
        "E.  The regression of analysis 4: the nano map predicted from receptor "
        "mRNA and synaptic density",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    return saved_working(fig, save)


def draw_map(
    fig: plt.Figure,
    ax: plt.Axes,
    img: np.ndarray,
    cmap: str,
    span: float | None,
    floor: float,
    n_structures: int,
    labels: tuple[str | None, str | None, bool],
) -> None:
    """Draw one map of F_maps, with its column's colour bar under the last plane.

    `span` is the symmetric limit of the leftover, None for a map of ranks, whose
    scale starts `floor` of the ranks below rank 1 so that nothing draws as black.
    `labels` holds the column's title (first plane), the plane's label (first
    column), and whether the colour bar goes under it (last plane).
    """
    title, ylabel, with_bar = labels
    if span is None:
        lo = 1 - floor * (n_structures - 1)
        im = ax.imshow(
            img, cmap=cmap, vmin=lo, vmax=n_structures, interpolation="nearest"
        )
    else:
        im = ax.imshow(img, cmap=cmap, vmin=-span, vmax=span, interpolation="nearest")

    # the image rasterised in the EPS, the text kept vector; NaN is
    # transparent, so the black face is the ground
    im.set_rasterized(True)
    ax.set_facecolor("black")
    ax.set_xticks([])
    ax.set_yticks([])
    for side in ax.spines.values():
        side.set_visible(False)
    if title is not None:
        ax.set_title(title, fontsize=9)
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=8)
    if with_bar:
        cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02, orientation="horizontal")
        cb.ax.tick_params(labelsize=6.5)
        if span is None:
            label = "rank among structures"
        else:
            label = "observed minus predicted (ranks)"
        cb.set_label(label, fontsize=7)
        if span is None:
            # the floor sits below rank 1 so that nothing draws as black;
            # the ticks still stop at the real range
            cb.set_ticks([1, 25, 50, 75, 100, n_structures])


def plot_regression_maps(
    planes: list[dict],
    observed: np.ndarray,
    predicted: np.ndarray,
    leftover: np.ndarray,
    structures: list[str],
    floor: float,
    save: Path | None = None,
) -> plt.Figure:
    """F_maps: observed, predicted and residual maps on coronal planes.

    `planes` holds per plane its parcellation indices (lab), its number in the full
    CCF at 10 um (ccf_plane) and the structure names of the indices (names).
    """
    # (title, value per structure, colormap, symmetric limit or None for ranks)
    maps = [
        ("the nano map", dict(zip(structures, observed)), "hot", None),
        (
            "predicted from Gria1\nand synapse density",
            dict(zip(structures, predicted)),
            "hot",
            None,
        ),
        (
            "what is left over",
            dict(zip(structures, leftover)),
            "RdBu_r",
            float(np.abs(leftover).max()),
        ),
    ]
    fig, axes = plt.subplots(len(planes), 3, figsize=(10.2, 3.2 * len(planes)))
    axes = np.atleast_2d(axes)
    for r, plane in enumerate(planes):
        # the plane in the full CCF at 10 um, as the route's other figures give it
        ylabel = f"CCF plane {plane['ccf_plane']} / 10 um"
        for c, (title, values, cmap, span) in enumerate(maps):
            img = paint(plane["lab"], values, plane["names"])
            labels = (
                title if r == 0 else None,
                ylabel if c == 0 else None,
                r == len(planes) - 1,
            )
            draw_map(fig, axes[r, c], img, cmap, span, floor, len(structures), labels)
    fig.suptitle(
        "F.  The same three quantities on the brain.  Red in the third column is "
        "a higher nano rank than\nGria1 and synapse density predict, blue is lower.  "
        "Black is outside the brain, or a structure the\nfit does not use.",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return saved_working(fig, save)
