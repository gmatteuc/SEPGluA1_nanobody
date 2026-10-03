"""The summary figures of the beyond-abundance result, with their statistics.

adult.beyond_density and adult.beyond_controls draw working figures, seven of them,
for checking each step. This module draws the version to show: four panels in the
palette of the young-against-adult figures, every number with an interval, and an
EPS beside every PNG for editing.

    A  what explains the map         cross-validated R2 of each explanation against
                                     the ceiling, with bootstrap intervals
    B  the leftover is real          the two half-cohort agreements against the
                                     null of "the leftover is noise", with a p
    C  where the leftover lives      the structures that carry it
    D  seven ways it could be wrong  the controls, at a glance

Intervals come from resampling structures, not animals. The claim is about where in
the brain the map departs from prediction, so what would differ in a repeat of the
analysis is which structures were measurable. Each bootstrap replicate draws 125
structures with replacement and recomputes the whole quantity; the uncertainty over
animals is carried by the ceiling and the half-cohort split, a separate axis shown
separately.

The p of panel B asks whether the leftover could be noise. If it were, the leftover
of one half-cohort would say nothing about that of the other, and their correlation
would sit at zero. The null shuffles which structure is which in one half, which
breaks the correspondence while keeping both distributions exactly as they are.
R2 is cross-validated throughout, because a model with more terms always fits
better on the data it was fitted to.

Writes, in adult_v2/beyond/for_sami/ under the data root:

    A_what_explains.png/.eps      B_leftover_real.png/.eps
    C_where.png/.eps              D_controls.png/.eps
    numbers_for_the_caption.txt   every figure's numbers as a sentence

Run by run_beyond_figures.py.
"""

import csv
from collections.abc import Sequence

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import to_rgb
from scipy.stats import rankdata, spearmanr

from sepmap.adult.beyond_density import (
    ADULTS,
    OUT,
    build_covariates,
    cv_r2,
    flexible,
    half_map,
    half_splits,
    prepare,
    residual,
    spearman_brown,
)
from sepmap.config import SETTINGS
from sepmap.plotting import DARK_GREY, MID_GREY, RED, tidy

# the bootstrap replicates, the permutations of the noise null, and the half-cohort
# splits used inside each replicate
BEYOND_FIGURES = SETTINGS["beyond_figures"]

FIGS = OUT / "for_sami"

# the palette of the young-against-adult figures, so the two sets read as one: red
# for what is shown, dark grey for its comparison, mid grey for context (plotting),
# and this blue only for the signed panel, where red-blue has a zero
BLUE = "#2e5f8a"


def save(fig: plt.Figure, name: str) -> None:
    """Save `fig` in FIGS as a PNG at 220 dpi and an EPS to edit, then close it."""
    FIGS.mkdir(parents=True, exist_ok=True)
    png = FIGS / (name + ".png")
    fig.savefig(png, dpi=220)
    fig.savefig(FIGS / (name + ".eps"), format="eps")
    plt.close(fig)
    print(f"  -> {png}  (+ .eps)")


def percentile_interval(values: Sequence[float] | np.ndarray) -> tuple[float, float]:
    """The 95% interval of `values`, from the 2.5th to the 97.5th percentile."""
    return float(np.percentile(values, 2.5)), float(np.percentile(values, 97.5))


def bootstrap_models(
    y: np.ndarray,
    models: list[tuple[str, list[np.ndarray]]],
    ceiling: float,
    rng: np.random.Generator,
) -> tuple[list[float], list[tuple[float, float]]]:
    """Cross-validated R2 of each model, with an interval, resampling structures.

    Each of the beyond_figures.n_boot replicates ranks `y` and every predictor
    again over the
    structures it drew. Returns the point values and the (low, high) intervals.
    """
    point = [cv_r2(y, xs) for _, xs in models]
    n = len(y)
    draws = np.empty((BEYOND_FIGURES["n_boot"], len(models)))
    for b in range(BEYOND_FIGURES["n_boot"]):
        take = rng.integers(0, n, n)
        yb = rankdata(y[take])
        for j, (_, xs) in enumerate(models):
            draws[b, j] = cv_r2(yb, [rankdata(x[take]) for x in xs])
    return point, [percentile_interval(draws[:, j]) for j in range(len(models))]


def bootstrap_replication(
    nano: dict[str, dict[str, float]],
    structures: list[str],
    predictors: list[np.ndarray],
    splits: list[tuple[list[int], list[int]]],
    rng: np.random.Generator,
) -> tuple[float, tuple[float, float]]:
    """Half-against-half agreement of the leftover, with an interval.

    The point value uses every split; the beyond_figures.n_boot_halves replicates
    each resample the structures and use the same boot_splits random splits.
    """

    def agreement(idx, use_splits):
        """Mean agreement of the two leftovers over `use_splits`, on structures `idx`."""
        picked = [structures[i] for i in idx]
        return float(
            np.mean(
                [
                    spearmanr(
                        residual(half_map(nano, a, picked), [p[idx] for p in predictors]),
                        residual(half_map(nano, b, picked), [p[idx] for p in predictors]),
                    ).statistic
                    for a, b in use_splits
                ]
            )
        )

    full = np.arange(len(structures))
    point = agreement(full, splits)
    few = [
        splits[i]
        for i in rng.choice(len(splits), BEYOND_FIGURES["boot_splits"], replace=False)
    ]
    draws = [
        agreement(rng.integers(0, len(structures), len(structures)), few)
        for _ in range(BEYOND_FIGURES["n_boot_halves"])
    ]
    return point, percentile_interval(draws)


def noise_null(
    nano: dict[str, dict[str, float]],
    structures: list[str],
    predictors: list[np.ndarray],
    splits: list[tuple[list[int], list[int]]],
    rng: np.random.Generator,
) -> tuple[float, np.ndarray, float]:
    """How well the two halves would agree if the leftover were noise, and the p.

    Shuffling which structure is which in one half breaks the correspondence
    between the two halves while leaving both leftovers exactly as they are, so
    the null says "these two are unrelated" and nothing else. Uses the first
    split; returns the observed agreement, the beyond_figures.n_perm null values and
    the
    two-sided p.
    """
    a, b = splits[0]
    ra = residual(half_map(nano, a, structures), predictors)
    rb = residual(half_map(nano, b, structures), predictors)
    observed = float(spearmanr(ra, rb).statistic)
    null = np.empty(BEYOND_FIGURES["n_perm"])
    for i in range(BEYOND_FIGURES["n_perm"]):
        null[i] = spearmanr(ra, rng.permutation(rb)).statistic
    p = float(
        (np.sum(np.abs(null) >= abs(observed)) + 1) / (BEYOND_FIGURES["n_perm"] + 1)
    )
    return observed, null, p


# ===== The panels =====


def panel_a(
    point: list[float],
    intervals: list[tuple[float, float]],
    labels: list[str],
    ceiling: float,
    ceiling_ci: tuple[float, float],
) -> None:
    """Draw panel A: each explanation's cross-validated R2 against the ceiling."""
    fig, ax = plt.subplots(figsize=(7.8, 4.4))
    y = np.arange(len(labels))
    lo = [p - i[0] for p, i in zip(point, intervals)]
    hi = [i[1] - p for p, i in zip(point, intervals)]
    colours = [MID_GREY] * (len(labels) - 1) + [RED]
    ax.barh(y, point, color=colours, edgecolor="0.25", linewidth=0.5)
    ax.errorbar(
        point, y, xerr=[lo, hi], fmt="none", ecolor="0.2", elinewidth=0.9, capsize=2.5
    )
    ax.axvline(ceiling, color=DARK_GREY, lw=1.8)

    # the ceiling's interval in the grey that 15% of DARK_GREY gives on white, opaque
    # and under the bars: EPS has no transparency and printed the band solid
    band = 0.85 * np.ones(3) + 0.15 * np.array(to_rgb(DARK_GREY))
    ax.axvspan(ceiling_ci[0], ceiling_ci[1], color=band, lw=0, zorder=0)
    ax.text(
        ceiling - 0.015,
        0.97,
        "ceiling: all of the map\nthat is explainable",
        transform=ax.get_xaxis_transform(),
        color=DARK_GREY,
        fontsize=8,
        ha="right",
        va="top",
    )
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8.5)
    ax.set_xlim(0, 1.02)
    ax.set_xlabel(
        "variance of the surface-GluA1 map predicted\n"
        "(cross-validated, 95% bootstrap interval over structures)",
        fontsize=8.5,
    )
    ax.set_title(
        "A.  Receptor abundance and synaptic density explain part of the map",
        fontsize=10,
        loc="left",
    )
    tidy(ax)
    fig.tight_layout()
    save(fig, "A_what_explains")


def panel_b(
    map_agreement: list[float],
    leftover_agreement: list[float],
    null: np.ndarray,
    p: float,
    rep_point: float,
    rep_ci: tuple[float, float],
) -> None:
    """Draw panel B: the agreement of the map, of the leftover and of the noise null."""
    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    bins = np.linspace(-0.45, 1.0, 120)
    for values, colour, label in (
        (null, MID_GREY, "if the leftover were noise"),
        (map_agreement, DARK_GREY, "the map itself"),
        (leftover_agreement, RED, "what is left of it"),
    ):
        counts, edges = np.histogram(values, bins=bins)
        ax.bar(
            edges[:-1],
            counts / counts.max(),
            width=np.diff(edges),
            align="edge",
            color=colour,
            edgecolor="none",
            label=label,
        )
    ax.set_xlabel(
        "agreement between two independent halves of the cohort (Spearman)",
        fontsize=8.5,
    )
    ax.set_ylabel("how often, each curve scaled to its own peak", fontsize=8.5)
    ax.set_ylim(0, 1.15)
    ax.legend(fontsize=8, frameon=False, loc="upper left")
    ax.set_title(
        f"B.  What is left over replicates across animals\n"
        f"leftover {rep_point:.3f} [{rep_ci[0]:.3f}, {rep_ci[1]:.3f}], "
        f"against noise p < {max(p, 1e-4):.0e}",
        fontsize=10,
        loc="left",
    )
    tidy(ax)
    fig.tight_layout()
    save(fig, "B_leftover_real")


def panel_c(res: np.ndarray, structures: list[str], n_show: int = 9) -> None:
    """Draw panel C: the `n_show` largest residuals of each sign."""
    order = np.argsort(-res)
    show = list(order[:n_show]) + list(order[-n_show:])
    fig, ax = plt.subplots(figsize=(7.6, 5.6))
    pos = np.arange(len(show))
    ax.barh(
        pos,
        [res[i] for i in show],
        color=[RED if res[i] > 0 else BLUE for i in show],
        edgecolor="0.25",
        linewidth=0.4,
    )
    ax.set_yticks(pos)
    ax.set_yticklabels([structures[i][:42] for i in show], fontsize=8)
    ax.invert_yaxis()
    ax.axvline(0, color="0.3", lw=0.8)
    ax.set_xlabel(
        "surface GluA1, minus what abundance and density predict (ranks)", fontsize=8.5
    )
    ax.set_title("C.  The leftover is anatomically organised", fontsize=10, loc="left")
    tidy(ax)
    fig.tight_layout()
    save(fig, "C_where")


def panel_d(controls: list[dict[str, str]]) -> None:
    """Draw panel D: the rows of controls.csv as a table of verdicts."""
    fig, ax = plt.subplots(figsize=(9.4, 4.0))
    ax.axis("off")
    ax.set_title(
        "D.  Seven ways the result could be an artefact, and the number for each",
        fontsize=10,
        loc="left",
    )
    for i, row in enumerate(controls):
        y = 1 - (i + 1) / (len(controls) + 1)
        ok = row["verdict"] == "pass"
        ax.text(
            0.0, y, row["control"], fontsize=9, va="center", color="0.15" if ok else RED
        )
        ax.text(0.38, y, row["number"][:78], fontsize=8, va="center", color="0.35")
        ax.text(
            0.985,
            y,
            "ruled out" if ok else "CHECK",
            fontsize=9,
            va="center",
            ha="right",
            color=DARK_GREY if ok else RED,
        )
    fig.tight_layout()
    save(fig, "D_controls")


def ceiling_with_interval(
    nano: dict[str, dict[str, float]],
    structures: list[str],
    splits: list[tuple[list[int], list[int]]],
    rng: np.random.Generator,
) -> tuple[list[float], float, tuple[float, float]]:
    """The ceiling, with an interval over structures.

    Returns the half-cohort agreement of every split, the ceiling and its interval
    from beyond_figures.n_boot_halves replicates, each drawn from `rng`.
    """
    map_agreement = [
        spearmanr(half_map(nano, a, structures), half_map(nano, b, structures)).statistic
        for a, b in splits
    ]
    ceiling = spearman_brown(float(np.mean(map_agreement))) ** 2
    boot_ceiling = []
    for _ in range(BEYOND_FIGURES["n_boot_halves"]):
        take = rng.integers(0, len(structures), len(structures))
        picked = [structures[i] for i in take]
        few = [
            splits[i]
            for i in rng.choice(len(splits), BEYOND_FIGURES["boot_splits"], replace=False)
        ]
        boot_ceiling.append(
            spearman_brown(
                float(
                    np.mean(
                        [
                            spearmanr(
                                half_map(nano, a, picked), half_map(nano, b, picked)
                            ).statistic
                            for a, b in few
                        ]
                    )
                )
            )
            ** 2
        )
    ceiling_ci = percentile_interval(boot_ceiling)
    print(f"ceiling {ceiling:.3f} [{ceiling_ci[0]:.3f}, {ceiling_ci[1]:.3f}]")
    return map_agreement, ceiling, ceiling_ci


def explanations(
    covariates: dict[str, np.ndarray],
    model: list[np.ndarray],
    y: np.ndarray,
    ceiling: float,
    rng: np.random.Generator,
) -> tuple[list[tuple[str, list[np.ndarray]]], list[float], list[tuple[float, float]]]:
    """What each explanation predicts, with intervals; printed.

    Returns the models as (label, predictors), their cross-validated R2 and
    intervals.
    """
    c = covariates
    models = [
        ("receptor abundance\n(Gria1-4)", [c["abundance"]]),
        ("synaptic markers\n(11 pre- and postsynaptic)", [c["markers"]]),
        ("postsynaptic gene\nexpression (188 genes)", [c["psd_pc1"]]),
        ("tissue autofluorescence\n(the same brains)", [c["autofluo"]]),
        ("all four together", model),
    ]
    point, intervals = bootstrap_models(y, models, ceiling, rng)
    for (label, _), pt, iv in zip(models, point, intervals):
        print(
            f"  {label.splitlines()[0]:34s} CV R2 {pt:.3f} [{iv[0]:.3f}, {iv[1]:.3f}]"
            f"   {pt / ceiling:5.1%} of the ceiling"
        )
    return models, point, intervals


def leftover_replication(
    nano: dict[str, dict[str, float]],
    structures: list[str],
    model: list[np.ndarray],
    splits: list[tuple[list[int], list[int]]],
    rng: np.random.Generator,
) -> tuple[list[float], float, tuple[float, float], np.ndarray, float]:
    """The leftover's replication, its interval and the noise null; printed.

    Returns the agreement of every split, the replication and its interval, the
    null values and the p.
    """
    leftover_agreement = [
        spearmanr(
            residual(half_map(nano, a, structures), model),
            residual(half_map(nano, b, structures), model),
        ).statistic
        for a, b in splits
    ]
    rep_point, rep_ci = bootstrap_replication(nano, structures, model, splits, rng)
    observed, null, p = noise_null(nano, structures, model, splits, rng)
    print(
        f"  leftover replicates {rep_point:.3f} [{rep_ci[0]:.3f}, {rep_ci[1]:.3f}]; "
        f"against the noise null p = {p:.2e}"
    )
    return leftover_agreement, rep_point, rep_ci, null, p


def load_controls() -> list[dict[str, str]]:
    """The rows of run_beyond_controls' controls.csv, which panel D needs, if any."""
    controls = []
    path = OUT / "controls.csv"
    if path.exists():
        with open(path, newline="", encoding="utf-8") as fh:
            controls = list(csv.DictReader(fh))
    return controls


def write_caption_numbers(
    ceiling: float,
    ceiling_ci: tuple[float, float],
    point: list[float],
    unexplained: float,
    splits: list[tuple[list[int], list[int]]],
    map_agreement: list[float],
    rep_point: float,
    rep_ci: tuple[float, float],
    p: float,
) -> None:
    """Write numbers_for_the_caption.txt: every figure's numbers as a sentence."""
    lines = [
        "Numbers for the captions (all on 125 grey-matter structures, ten adult mice).",
        "",
        f"A. The map is reproducible enough that {ceiling:.1%} of its variance is",
        f"   explainable in principle (95% CI {ceiling_ci[0]:.1%} "
        f"to {ceiling_ci[1]:.1%}).",
        f"   Receptor abundance predicts {point[0] / ceiling:.0%} "
        f"of that, synaptic markers",
        f"   {point[1] / ceiling:.0%}, postsynaptic gene expression "
        f"{point[2] / ceiling:.0%},",
        f"   tissue autofluorescence {max(point[3], 0) / ceiling:.0%}. "
        f"All four together,",
        f"   allowed a non-linear relationship, predict {point[-1] / ceiling:.0%},",
        f"   leaving {unexplained:.0%} unexplained.",
        "",
        f"B. Splitting the ten animals into two fives every possible way "
        f"({len(splits)} splits),",
        f"   the two half-cohort maps agree at {np.mean(map_agreement):.3f}. "
        f"After the four",
        f"   explanations are removed, the two leftovers still agree at {rep_point:.3f}",
        f"   (95% CI {rep_ci[0]:.3f} to {rep_ci[1]:.3f}). Were the leftover noise, that",
        f"   agreement would be zero; p < {max(p, 1e-4):.0e} by permutation.",
        "",
        "C. Structures where surface GluA1 most exceeds and falls short of what",
        "   abundance and density predict, in ranks among the 125.",
        "",
        "D. Seven controls, each ruling out a way the leftover could be an artefact.",
        "",
        "What this does NOT show: that the leftover is the surface fraction. A residual",
        "is only ever what the model left out.",
    ]
    with open(FIGS / "numbers_for_the_caption.txt", "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"  -> {FIGS / 'numbers_for_the_caption.txt'}")


def main() -> None:
    """Compute the statistics, draw panels A to D and write the caption numbers.

    One generator feeds every interval and the null, in this order: the ceiling,
    the explanations, the leftover's replication, the noise null.
    """
    FIGS.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(0)

    # the structures and the quoted model, as adult.beyond_density builds them
    nano, _, expr, role, auto, structures = prepare()
    covariates, _, _ = build_covariates(nano, expr, role, auto, structures)
    splits = half_splits()
    y = half_map(nano, range(len(ADULTS)), structures)
    model = flexible(list(covariates.values()))
    print(
        f"{len(structures)} structures, {len(ADULTS)} adults, {len(splits)} half-splits"
    )

    # the ceiling, what each explanation predicts, and how the leftover replicates
    map_agreement, ceiling, ceiling_ci = ceiling_with_interval(
        nano, structures, splits, rng
    )
    models, point, intervals = explanations(covariates, model, y, ceiling, rng)
    unexplained = 1 - point[-1] / ceiling
    print(f"  unexplained share of the ceiling: {unexplained:.1%}")
    leftover_agreement, rep_point, rep_ci, null, p = leftover_replication(
        nano, structures, model, splits, rng
    )

    res = residual(y, model)

    # the panels; panel D needs run_beyond_controls' table
    controls = load_controls()
    panel_a(point, intervals, [m[0] for m in models], ceiling, ceiling_ci)
    panel_b(map_agreement, leftover_agreement, null, p, rep_point, rep_ci)
    panel_c(res, structures)
    if controls:
        panel_d(controls)

    # the numbers for the captions
    write_caption_numbers(
        ceiling,
        ceiling_ci,
        point,
        unexplained,
        splits,
        map_agreement,
        rep_point,
        rep_ci,
        p,
    )
