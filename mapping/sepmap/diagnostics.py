"""Diagnostic sheets for every step of the route, so it can be audited, not believed.

One figure per question, written to comparisons_v2/processing_diagnostics/ under
the data root, with an index README beside them:

    01_tissue_<mouse>.png       what was counted as tissue, drawn on the brain
    02_levels_<mouse>.png       where the background and the threshold sit in the
                                intensity distributions
    03_coverage.png             how much of the atlas each brain covers, plane by
                                plane: where a cohort mean rests on few mice
    04_warp_<mouse>.png         the same plane before and after DeMBA -> CCF, and
                                every region's value before against after
    05_cohort_n.png             how many brains contribute at each voxel
    06_scaling.png              the per-mouse scale factors, and what they do to
                                the cortex distributions
    07_route_agreement.png      voxelwise (warped) against region-wise (never
                                warped), structure by structure
    08_mask_vs_p6bis.png        this route's tissue mask against the one of
                                run_normalise_groups / run_nano_equalisation,
                                checked by eye, where the latter works
    09_denominators.png         autofluorescence and SEP as reference channels,
                                and whether a structure's answer depends on which
                                one is used

The per-brain sheets run_add_sep_channel writes when it carries the SEP channel
into registered space sit beside these, in processing_diagnostics/sep_channel/.

Run by run_diagnostics.py, after a full pass, or with mouse names to refresh a
few sheets.
"""

import csv

import h5py
import matplotlib.pyplot as plt
import numpy as np

from sepmap.config import SETTINGS
from sepmap.plotting import DARK_GREY, RED, save_figure
from sepmap.volumes.cohort import COHORTS, PER_MOUSE_CCF
from sepmap.volumes.cohort import OUT_ROOT as CCF_ROOT
from sepmap.volumes.per_mouse import (
    DATA,
    MICE,
    annotation_20,
    isocortex_ids,
)
from sepmap.volumes.per_mouse import OUT as PER_MOUSE

# the tissue mask's threshold, in MADs above the auto background
TISSUE = SETTINGS["tissue"]

OUT = DATA / "comparisons_v2" / "processing_diagnostics"

# the young brains are on their own DeMBA atlas, the adults on the CCF
YOUNG = [m for m, v in MICE.items() if v[1] != "ccf"]
ADULT = [m for m, v in MICE.items() if v[1] == "ccf"]


def show(
    ax: plt.Axes, img: np.ndarray, mask: np.ndarray | None = None, p: float = 99.5
) -> None:
    """A plane, dorsal up and ventral down, scaled to its own tissue.

    A plane of these volumes is (DV, ML), so it is drawn as it comes: rows run
    dorsal to ventral, columns left to right. The slice figures elsewhere use
    the same convention.
    """
    # grey up to the p-th percentile of the plane's positive values (1.0 if it has
    # none), so each plane fills the grey range whatever its brightness
    im = np.asarray(img, float)
    hi = np.percentile(im[im > 0], p) if (im > 0).any() else 1.0
    ax.imshow(np.clip(im / hi, 0, 1), cmap="gray", origin="upper")

    # the mask, if given, as a red outline
    if mask is not None:
        ax.contour(
            np.asarray(mask, float), levels=[0.5], colors="#e74c3c", linewidths=0.9
        )
    ax.axis("off")


def sheet_tissue(mouse: str, ann: np.ndarray, z: np.lib.npyio.NpzFile) -> None:
    """Sheet 01: the tissue mask on four planes, with the atlas outline for scale.

    `ann` is the brain's own atlas at 20 um and `z` its per-mouse file.
    """
    # the background-subtracted nano, the tissue mask and the atlas brain
    sig, tissue = z["sig"].astype(np.float32), z["tissue"]
    brain = ann > 0

    # four planes between the 5th and 95th percentiles of the planes where tissue
    # covers more than 20% of the atlas brain (no coverage, NaN, without atlas brain)
    cov = np.array(
        [
            tissue[k][brain[k]].mean() if brain[k].any() else np.nan
            for k in range(tissue.shape[0])
        ]
    )
    planes = (
        np.linspace(*np.percentile(np.flatnonzero(cov > 0.2), [5, 95]), 4)
        .round()
        .astype(int)
    )
    fig, axes = plt.subplots(1, 4, figsize=(19, 5.2))
    for ax, k in zip(axes, planes):
        # the plane with the mask in red and the atlas brain outlined in blue
        show(ax, sig[k], tissue[k])
        ax.contour(
            (ann[k] > 0).astype(float), levels=[0.5], colors="#3498db", linewidths=0.6
        )

        # the share of the atlas brain the mask covers on this plane
        ax.set_title(
            f"plane {k}   tissue {100 * tissue[k][brain[k]].mean():.0f}% of atlas brain",
            fontsize=10,
        )

    # the mask's rule and the two backgrounds subtracted
    fig.suptitle(
        f"{mouse}: red = tissue mask (auto channel above background "
        f"+ {TISSUE['mad_k']:g} MAD, "
        "and nano non-zero over half the voxel), "
        "blue = atlas brain.  Background subtracted: "
        f"nano {float(z['bg_nano']):.0f}, auto {float(z['bg_auto']):.0f} counts",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))

    # PNG only, as for every sheet: the EPS beside a figure is for a paper
    save_figure(fig, OUT / f"01_tissue_{mouse}.png", dpi=95, eps=False)
    plt.close(fig)


def draw_levels(
    ax: plt.Axes,
    arr: np.ndarray,
    name: str,
    bg: float,
    off: np.ndarray,
    tissue: np.ndarray,
    mad: float,
) -> None:
    """Draw one raw channel, tissue against off tissue, with its background marked.

    The auto channel also gets the mask threshold.
    """
    # off tissue in grey and tissue in red, the background added back to give raw
    # counts, up to four times the background; a subsample is plenty for a histogram
    ax.hist(
        arr[off][::17] + bg,
        bins=200,
        range=(0, 4 * bg),
        color="#95a5a6",
        label="off tissue",
        density=True,
    )
    ax.hist(
        arr[tissue][::37] + bg,
        bins=200,
        range=(0, 4 * bg),
        color=RED,
        alpha=0.6,
        label="tissue",
        density=True,
    )

    # the background, and for auto the mask threshold, tissue.mad_k MADs above it
    ax.axvline(bg, color="k", lw=1.2, label=f"background {bg:.0f}")
    if name == "auto":
        threshold = bg + TISSUE["mad_k"] * mad
        ax.axvline(
            threshold,
            color="#2980b9",
            lw=1.2,
            ls="--",
            label=f"mask threshold {threshold:.0f}",
        )
    ax.set_xlabel(f"{name} channel, raw counts")
    ax.set_ylabel("density")
    ax.legend(fontsize=8)
    ax.set_title(f"{name}: tissue vs off tissue", fontsize=10)


def draw_cortex_scaling(
    ax: plt.Axes,
    sig: np.ndarray,
    tissue: np.ndarray,
    ann: np.ndarray,
    z: np.lib.npyio.NpzFile,
) -> None:
    """Draw the isocortex after scaling by its mean, which should centre on 1."""
    # the isocortex tissue voxels over the brain's isocortex mean, every 37th
    iso = tissue & np.isin(ann, ISO)
    ax.hist(
        (sig[iso] / float(z["cortex_mean"]))[::37],
        bins=200,
        range=(0, 3),
        color=RED,
        density=True,
    )

    # a line at 1, where the scaled cortex should centre
    ax.axvline(1, color="k", lw=1.2)
    ax.set_xlabel("sig / isocortex mean")
    ax.set_ylabel("density")
    ax.set_title(
        f"cortex after scaling (mean {float(z['cortex_mean']):.0f} counts -> 1.0)",
        fontsize=10,
    )


def sheet_levels(mouse: str, ann: np.ndarray, z: np.lib.npyio.NpzFile) -> None:
    """Sheet 02: the intensity distributions the mask and the background rest on.

    `ann` is the brain's own atlas at 20 um and `z` its per-mouse file; the
    histograms take every 17th off-tissue and every 37th tissue voxel.
    """
    # the background-subtracted channels, the mask, and the scalars run_per_mouse
    # measured: the two backgrounds and the auto MAD
    sig, auto, tissue = (
        z["sig"].astype(np.float32),
        z["auto"].astype(np.float32),
        z["tissue"],
    )
    brain = ann > 0
    bg_n, bg_a, mad = float(z["bg_nano"]), float(z["bg_auto"]), float(z["mad_auto"])

    # imaged (raw auto above zero), outside the atlas brain
    off = (~brain) & (auto > -bg_a)

    # the raw auto and nano channels, tissue against off tissue
    fig, axes = plt.subplots(1, 3, figsize=(17, 4.6))
    for ax, (arr, name, bg) in zip(axes[:2], ((auto, "auto", bg_a), (sig, "nano", bg_n))):
        draw_levels(ax, arr, name, bg, off, tissue, mad)

    # the isocortex after scaling by its mean, which should centre on 1
    draw_cortex_scaling(axes[2], sig, tissue, ann, z)

    fig.suptitle(
        f"{mouse}: what the background subtraction and the mask threshold "
        "actually separate",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    save_figure(fig, OUT / f"02_levels_{mouse}.png", dpi=95, eps=False)
    plt.close(fig)


def line_colours(n: int) -> np.ndarray:
    """`n` distinguishable line colours with no green in them (plasma, trimmed)."""
    return plt.get_cmap("plasma")(np.linspace(0.0, 0.88, max(n, 2)))


def sheet_coverage() -> None:
    """Sheet 03: per brain, the fraction of the atlas brain with tissue, per plane."""
    # the young above, the adults below, each brain on the atlas it was measured on
    fig, axes = plt.subplots(2, 1, figsize=(14, 8))
    for ax, group, label in (
        (axes[0], YOUNG, "young, on each brain's own atlas"),
        (axes[1], ADULT, "adults, on the CCF"),
    ):
        cols = line_colours(len(group))
        for ci, mouse in enumerate(group):
            # the fraction of the atlas brain with tissue, plane by plane (NaN where
            # the atlas has no brain)
            ann = ANN[MICE[mouse][1]]
            brain = ann > 0
            t = np.load(PER_MOUSE / (mouse + ".npz"))["tissue"]
            cov = np.array(
                [
                    t[k][brain[k]].mean() if brain[k].any() else np.nan
                    for k in range(t.shape[0])
                ]
            )

            # a line per brain, its legend counting the planes over 20% covered
            ax.plot(
                np.arange(len(cov)),
                100 * cov,
                lw=1.3,
                color=cols[ci],
                label=f"{mouse} ({np.nansum(cov > 0.2):.0f} planes)",
            )
        ax.set_xlabel("atlas plane within the crop")
        ax.set_ylabel("% of atlas brain with tissue")
        ax.set_title(label, fontsize=10)
        ax.grid(lw=0.3, alpha=0.6)
        ax.legend(fontsize=7, ncol=2)

    fig.suptitle(
        "Coverage: where a cohort mean rests on every brain, "
        "and where it rests on one or two",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    save_figure(fig, OUT / "03_coverage.png", dpi=95, eps=False)
    plt.close(fig)


def sheet_warp(mouse: str) -> None:
    """Sheet 04: one young brain before and after the transform, drawn and in numbers.

    Two planes in the brain's own atlas beside the same planes in the CCF, and
    the mean of every third region before against after, for the regions with
    more than 300 tissue voxels on both sides.
    """
    # the brain on its own atlas (n), and carried into the CCF (c)
    ann_n = ANN[MICE[mouse][1]]
    zn = np.load(PER_MOUSE / (mouse + ".npz"))
    zc = np.load(PER_MOUSE_CCF / (mouse + ".npz"))
    sig_n, t_n = zn["sig"].astype(np.float32), zn["tissue"]
    sig_c, t_c = zc["sig"].astype(np.float32), zc["tissue"]

    # the same anatomy: match by the fraction of the covered range (planes with more
    # than 2000 tissue voxels), 30% and 60% of the way through
    cov_n = np.flatnonzero([t_n[k].sum() > 2000 for k in range(t_n.shape[0])])
    cov_c = np.flatnonzero([t_c[k].sum() > 2000 for k in range(t_c.shape[0])])
    fig = plt.figure(figsize=(17, 5.0))
    for i, f in enumerate((0.3, 0.6)):
        kn = int(cov_n[0] + f * (cov_n[-1] - cov_n[0]))
        kc = int(cov_c[0] + f * (cov_c[-1] - cov_c[0]))
        ax = fig.add_subplot(1, 3, i + 1)

        # both planes are (DV, ML): stacking along ML puts them side by side
        show(
            ax,
            np.concatenate([sig_n[kn], sig_c[kc]], axis=1),
            np.concatenate([t_n[kn], t_c[kc]], axis=1),
        )
        ax.set_title(
            f"{int(100 * f)}% through the stack: own atlas (left) | in CCF (right)",
            fontsize=10,
        )

    # region means in the own atlas against the CCF; the CCF annotation is the
    # adult crop, planes 90 to 539 of the full 20 um grid
    ax = fig.add_subplot(1, 3, 3)
    ccf_full = np.zeros(sig_c.shape, ANN["ccf"].dtype)
    ccf_full[90:540] = ANN["ccf"]

    # every third label of the own atlas, 0 (outside the brain) skipped, where both
    # sides have more than 300 tissue voxels
    a, b = [], []
    for idx in np.unique(ANN[MICE[mouse][1]])[1:][::3]:
        m1 = t_n & (ann_n == idx)
        m2 = t_c & (ccf_full == idx)
        if m1.sum() > 300 and m2.sum() > 300:
            a.append(sig_n[m1].mean())
            b.append(sig_c[m2].mean())

    # after against before on log axes, with the identity line a perfect transform
    # would put every region on
    a, b = np.array(a), np.array(b)
    ax.loglog(a, b, "o", ms=3, color=RED, alpha=0.6)
    lim = [min(a.min(), b.min()) * 0.9, max(a.max(), b.max()) * 1.1]
    ax.plot(lim, lim, "k-", lw=0.8)
    ax.set_xlabel("region mean, own atlas")
    ax.set_ylabel("region mean, in CCF")

    # the median change in log2, 0 for a transform that moves no signal between regions
    ax.set_title(
        f"{len(a)} regions, median |log2 change| {np.median(np.abs(np.log2(b / a))):.3f}",
        fontsize=10,
    )
    fig.suptitle(
        f"{mouse}: what the DeMBA -> CCF transform does (adults are placed, not warped)",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    save_figure(fig, OUT / f"04_warp_{mouse}.png", dpi=95, eps=False)
    plt.close(fig)


def sheet_cohort_n() -> None:
    """Sheet 05: how many brains support each voxel, per cohort, on four CCF planes."""
    # the cohorts the comparison uses, a row each, on four planes of the 20 um CCF grid
    cohorts = ["young", "young_P20", "adult"]
    planes = [150, 250, 350, 450]
    fig, axes = plt.subplots(
        len(cohorts), len(planes), figsize=(4.2 * len(planes), 3.6 * len(cohorts))
    )
    for r, cohort in enumerate(cohorts):
        # brains with a value at each voxel, on a scale up to the cohort's size
        n = np.load(CCF_ROOT / cohort / "cref_n.npy")
        for c, k in enumerate(planes):
            ax = axes[r, c]
            h = ax.imshow(
                n[k],
                cmap="magma",
                vmin=0,
                vmax=len(COHORTS[cohort]),
                origin="upper",
                interpolation="nearest",
            )
            # the plane numbered on the 10 um CCF grid, twice its 20 um index
            ax.set_title(f"{cohort}  CCF plane {2 * k} / 10 um", fontsize=9)
            ax.axis("off")
            plt.colorbar(h, ax=ax, fraction=0.03, pad=0.01)

    fig.suptitle(
        "How many brains contribute at each voxel (the n map the comparison thresholds)",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    save_figure(fig, OUT / "05_cohort_n.png", dpi=95, eps=False)
    plt.close(fig)


def sheet_scaling() -> None:
    """Sheet 06: the per-mouse numbers every normalisation rests on."""
    # per brain: mouse, cohort, nano and auto backgrounds, isocortex nano and auto
    rows = []
    for mouse in list(MICE):
        z = np.load(PER_MOUSE / (mouse + ".npz"))
        ann = ANN[MICE[mouse][1]]
        iso = z["tissue"] & np.isin(ann, ISO)
        rows.append(
            dict(
                mouse=mouse,
                cohort=MICE[mouse][0],
                bg_nano=float(z["bg_nano"]),
                bg_auto=float(z["bg_auto"]),
                cortex_nano=float(z["cortex_mean"]),
                cortex_auto=float(z["auto"].astype(np.float32)[iso].mean()),
            )
        )

    # the young and the adult brains, and both in that order along the x axis
    young = [r for r in rows if r["cohort"].startswith("young")]
    adult = [r for r in rows if not r["cohort"].startswith("young")]
    order = young + adult

    # nano background and isocortex mean per brain, the young first
    fig, axes = plt.subplots(1, 3, figsize=(17, 4.8))
    for ax, (column, name) in zip(
        axes,
        (
            ("bg_nano", "off-tissue background, nano"),
            ("cortex_nano", "isocortex mean, nano (bg-subtracted)"),
        ),
    ):
        # a dot per brain, red for the young, grey for the adults
        for i, r in enumerate(order):
            colour = RED if r["cohort"].startswith("young") else DARK_GREY
            ax.plot(i, r[column], "o", color=colour)

        # the brains by mouse ID, a dotted line between the young and the adults
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels(
            [r["mouse"].split("_")[0] for r in order], rotation=70, fontsize=7
        )
        ax.axvline(len(young) - 0.5, color="k", lw=0.6, ls=":")
        ax.set_ylabel("raw counts")
        ax.set_title(name, fontsize=10)
        ax.grid(lw=0.3, alpha=0.6)

    # isocortex nano against auto across brains
    ax = axes[2]
    for grp, col, lbl in ((young, RED, "young"), (adult, DARK_GREY, "adult")):
        ax.plot(
            [r["cortex_nano"] for r in grp],
            [r["cortex_auto"] for r in grp],
            "o",
            color=col,
            label=lbl,
        )
    ax.set_xlabel("isocortex mean, nano")
    ax.set_ylabel("isocortex mean, auto")

    # their correlation over every brain, young and adult together
    corr = np.corrcoef(
        [row["cortex_nano"] for row in rows], [row["cortex_auto"] for row in rows]
    )[0, 1]
    ax.set_title(f"the two channels across brains, r = {corr:.2f}", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(lw=0.3, alpha=0.6)

    # the spread of the adults' isocortex means, highest over lowest
    cortex = [r["cortex_nano"] for r in adult]
    fig.suptitle(
        "Per-brain levels: what the background subtraction removes, "
        "and what the cortex scaling divides by\n"
        f"The {max(cortex) / min(cortex):.0f}x spread among adults is why no analysis "
        "uses raw counts",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    save_figure(fig, OUT / "06_scaling.png", dpi=95, eps=False)
    plt.close(fig)


def sheet_route_agreement() -> None:
    """Sheet 07: the warped voxelwise route against the never-warped region route.

    The cref difference per structure from run_compare's region_table.csv
    against run_region_plot's region_stats.csv, for structures of at least 500
    voxels at 20 um. Skipped, with a printed line, when either table is missing.
    """
    d = DATA / "comparisons_v2" / "young_vs_adult"

    # the two tables, which run_compare and run_region_plot write
    missing = [
        name
        for name in ("region_table.csv", "region_stats.csv")
        if not (d / name).exists()
    ]
    if missing:
        print(
            f"sheet 07 skipped: {', '.join(missing)} missing in {d} "
            "(run run_compare.py and run_region_plot.py first)",
            flush=True,
        )
        return

    # acronym -> (log2 difference, division, voxels), and acronym -> difference
    with open(d / "region_table.csv", encoding="utf-8") as fh:
        vox = {
            r["acronym"]: (float(r["log2_cref"]), r["division"], int(r["voxels_20um"]))
            for r in csv.DictReader(fh)
        }
    with open(d / "region_stats.csv", encoding="utf-8") as fh:
        reg = {
            r["acronym"]: float(r["diff_log2"])
            for r in csv.DictReader(fh)
            if r["reading"] == "cref"
        }

    # the structures in both tables with at least 500 compared voxels, the isocortex
    # marked
    keys = [a for a in vox if a in reg and vox[a][2] >= 500]
    x = np.array([vox[a][0] for a in keys])
    y = np.array([reg[a] for a in keys])
    iso = np.array([vox[a][1] == "Isocortex" for a in keys])

    # one dot per structure, the isocortex in red, against the identity line
    fig, ax = plt.subplots(figsize=(7.2, 7))
    ax.plot(x[~iso], y[~iso], "o", ms=4, color="#95a5a6", label=f"other ({(~iso).sum()})")
    ax.plot(x[iso], y[iso], "o", ms=5, color=RED, label=f"isocortex ({iso.sum()})")
    lim = [min(x.min(), y.min()) - 0.1, max(x.max(), y.max()) + 0.1]
    ax.plot(lim, lim, "k-", lw=0.8)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel("log2 young/adult -- voxelwise, every brain warped into CCF")
    ax.set_ylabel("log2 young/adult -- region-wise, no warping")

    # the agreement, correlation and median absolute difference, over every structure
    # and over the isocortex alone
    ax.set_title(
        f"r = {np.corrcoef(x, y)[0, 1]:.3f} over {len(keys)} structures, "
        f"median |difference| {np.median(np.abs(x - y)):.3f} log2\n"
        f"isocortex only: r = {np.corrcoef(x[iso], y[iso])[0, 1]:.3f}, "
        f"median {np.median(np.abs(x[iso] - y[iso])):.3f}",
        fontsize=10,
    )
    ax.legend(fontsize=9)
    ax.grid(lw=0.3, alpha=0.6)
    fig.tight_layout()
    save_figure(fig, OUT / "07_route_agreement.png", dpi=110, eps=False)
    plt.close(fig)


def old_mask_plane(mask_4d, idx: int, k: int, mine: np.ndarray) -> np.ndarray:
    """Plane `k` of the run_normalise_groups mask of brain `idx`, at 20 um.

    The fraction of tissue in each 2x2x2 block of the 10 um mask, (DV, ML) like the
    per-mouse volumes.
    """
    # h5py reads the MATLAB array with its axes reversed, (mouse, ML, DV, AP) at 10 um;
    # 0 in the background mask is tissue
    old = np.zeros((mine.shape[1], mine.shape[2]), np.float32)
    for j in range(mine.shape[2]):
        # the tissue fraction over the block's two ML and two AP planes, along DV
        q = (
            sum(
                (np.asarray(mask_4d[idx, 2 * j + dj, :, 2 * k + dk]) == 0).astype(
                    np.float32
                )
                for dj in (0, 1)
                for dk in (0, 1)
            )
            / 4
        )

        # then over each pair of DV voxels, to the 20 um column
        old[:, j] = q.reshape(mine.shape[1], 2).mean(axis=1)
    return old


def draw_mask_plane(
    ax: plt.Axes,
    mouse: str,
    k: int,
    sig: np.ndarray,
    mine: np.ndarray,
    old: np.ndarray,
    ann: np.ndarray,
) -> float:
    """Draw plane `k` with both masks and the atlas brain; return the masks' Dice."""
    # the nano plane with this route's mask in red
    show(ax, sig[k], mine[k])

    # the old mask dashed blue, more than half tissue, and the atlas brain grey; both
    # are (DV, ML) like the plane underneath, so neither is transposed
    ax.contour(
        (old > 0.5).astype(float),
        levels=[0.5],
        colors="#3498db",
        linewidths=0.9,
        linestyles="--",
    )
    ax.contour(
        (ann[k] > 0).astype(float),
        levels=[0.5],
        colors="#cccccc",
        linewidths=0.6,
    )

    # the Dice coefficient of the two masks on this plane, 0 when both are empty
    agree = 2 * (mine[k] & (old > 0.5)).sum() / max(mine[k].sum() + (old > 0.5).sum(), 1)
    ax.set_title(f"{mouse}  plane {k}   Dice {agree:.3f}", fontsize=10)
    return float(agree)


def sheet_mask_vs_p6bis() -> None:
    """Sheet 08: this route's mask against run_normalise_groups', on two brains.

    The run_normalise_groups mask was checked by eye over many sessions, so it
    is the right thing to be measured against. It is not used here for two
    reasons: it is computed on the nano channel, which is the quantity being
    compared and is four times dimmer in pups; and its centre-of-mass guard
    returns an all-background mask for any partial section. Where it does work
    the two masks agree to better than 2% of voxels. Its 10 um mask is brought
    to 20 um as the fraction of tissue in each 2x2x2 block.
    """
    # mouse, its group's background mask file, the mouse's index in it, planes
    cases = [
        (
            "MG903_SepGluA_P20",
            DATA / "young" / "nano_4d_normalized_bkgmask_P20.mat",
            1,
            (200, 260, 330),
        ),
        (
            "CGF027_Gria1",
            DATA / "naive" / "nano_4d_normalized_bkgmask.mat",
            0,
            (150, 230, 300),
        ),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(18, 10.5))
    dice = []
    for row, (mouse, mask_file, idx, planes) in enumerate(cases):
        # this route's mask and nano for the brain, a row of the sheet
        ann = ANN[MICE[mouse][1]]
        z = np.load(PER_MOUSE / (mouse + ".npz"))
        mine, sig = z["tissue"], z["sig"].astype(np.float32)

        # the old mask from the HDF5 (version 7.3) .mat file, read piece by piece
        # rather than loaded whole
        with h5py.File(mask_file, "r") as f:
            mask_4d = f["recomputed_bkg_mask_4d"]
            for col, k in enumerate(planes):
                old = old_mask_plane(mask_4d, idx, k, mine)
                ax = axes[row, col]
                dice.append(draw_mask_plane(ax, mouse, k, sig, mine, old, ann))

    # the range of the Dice over the planes shown
    fig.suptitle(
        "red = v2 mask (auto channel), "
        "blue dashed = run_normalise_groups mask (nano channel), grey = atlas brain.\n"
        f"The two masks agree at Dice {min(dice):.3f} to {max(dice):.3f} on the "
        f"{len(dice)} planes shown",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    save_figure(fig, OUT / "08_mask_vs_p6bis.png", dpi=95, eps=False)
    plt.close(fig)


def denominator_rows() -> list[dict]:
    """Per brain with a SEP channel: mouse, cohort, isocortex nano, auto and SEP."""
    rows = []
    for mouse in list(MICE):
        # only the brains whose per-mouse file carries a SEP channel
        z = np.load(PER_MOUSE / (mouse + ".npz"))
        if "sep" not in z.files:
            continue

        # each channel's mean over the isocortex tissue, nano's as run_per_mouse
        # stored it
        ann = ANN[MICE[mouse][1]]
        iso = z["tissue"] & np.isin(ann, ISO)
        rows.append(
            dict(
                mouse=mouse,
                cohort=MICE[mouse][0],
                cortex_nano=float(z["cortex_mean"]),
                cortex_auto=float(z["auto"].astype(np.float32)[iso].mean()),
                cortex_sep=float(z["sep"].astype(np.float32)[iso].mean()),
            )
        )
    return rows


def panel_denominators(ax: plt.Axes, order: list[dict], young: list[dict]) -> None:
    """Draw what each denominator does with age, in the cortex, per brain."""
    # per brain, auto as an open circle and SEP as a filled square, red for the
    # young; labelled on the first brain only, so the legend has one entry each
    for i, r in enumerate(order):
        col = RED if r["cohort"].startswith("young") else DARK_GREY
        ax.plot(
            i,
            r["cortex_auto"],
            "o",
            color=col,
            mfc="none",
            label="auto" if i == 0 else None,
        )
        ax.plot(i, r["cortex_sep"], "s", color=col, label="SEP" if i == 0 else None)

    # the brains by mouse ID, a dotted line between the young and the adults
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([r["mouse"].split("_")[0] for r in order], rotation=70, fontsize=7)
    ax.axvline(len(young) - 0.5, color="k", lw=0.6, ls=":")
    ax.set_ylabel("isocortex mean, raw counts")
    ax.legend(fontsize=8)
    ax.set_title("the two denominators (open = auto, filled = SEP)", fontsize=10)
    ax.grid(lw=0.3, alpha=0.6)


def panel_nano_sep(ax: plt.Axes, young: list[dict], adult: list[dict]) -> None:
    """Draw how much of the nano difference SEP would absorb, across brains."""
    # a dot per brain, the young in red and the adults in grey
    for grp, col, lbl in ((young, RED, "young"), (adult, DARK_GREY, "adult")):
        ax.plot(
            [r["cortex_nano"] for r in grp],
            [r["cortex_sep"] for r in grp],
            "o",
            color=col,
            label=lbl,
        )
    ax.set_xlabel("isocortex mean, nano")
    ax.set_ylabel("isocortex mean, SEP")
    ax.set_title("nano against SEP across brains", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(lw=0.3, alpha=0.6)


def panel_swap(ax: plt.Axes) -> None:
    """Draw whether the young-adult difference survives the swap of denominator."""
    # each structure's young-adult log2 difference with each denominator, from
    # run_region_plot's region_stats.csv when it exists
    stats = DATA / "comparisons_v2" / "young_vs_adult" / "region_stats.csv"
    pairs = {}
    if stats.exists():
        with open(stats, encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                if r["reading"] in ("ratio", "sepratio"):
                    pairs.setdefault(r["acronym"], {})[r["reading"]] = float(
                        r["diff_log2"]
                    )

    # the structures with both, on square axes symmetric about zero, at least
    # 0.1 log2 each way
    keys = [k for k, v in pairs.items() if len(v) == 2]
    if keys:
        x = np.array([pairs[k]["ratio"] for k in keys])
        y = np.array([pairs[k]["sepratio"] for k in keys])
        lim = float(max(np.abs(np.concatenate([x, y])).max(), 0.1)) * 1.1

        # the identity line and the zero axes light, under the points
        ax.plot([-lim, lim], [-lim, lim], "-", color="#bbbbbb", lw=1)
        ax.axhline(0, color="#dddddd", lw=0.8)
        ax.axvline(0, color="#dddddd", lw=0.8)
        ax.plot(x, y, "o", ms=4, color=RED)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_title(
            f"young - adult, {len(keys)} structures (r = {np.corrcoef(x, y)[0, 1]:.2f})",
            fontsize=10,
        )
    else:
        # without the table, the step to run instead
        ax.text(
            0.5,
            0.5,
            "run run_region_plot.py first",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        ax.set_title("young - adult per structure", fontsize=10)
    ax.set_xlabel("log2 difference, nano / auto")
    ax.set_ylabel("log2 difference, nano / SEP")
    ax.grid(lw=0.3, alpha=0.6)


def sheet_denominators() -> None:
    """Sheet 09: the two reference channels side by side, and whether the answer moves.

    Only the brains that carry a SEP channel; the right panel needs
    run_region_plot's region_stats.csv and says so when it is missing.
    """
    # per brain: mouse, cohort, isocortex nano, auto and SEP
    rows = denominator_rows()
    if not rows:
        print("09 skipped: no brain carries a SEP channel yet", flush=True)
        return
    young = [r for r in rows if r["cohort"].startswith("young")]
    adult = [r for r in rows if not r["cohort"].startswith("young")]
    order = young + adult

    fig, axes = plt.subplots(1, 3, figsize=(17, 4.8))

    # what each denominator does with age; how much of the nano difference each
    # one would absorb; and the part that matters, does the difference survive
    panel_denominators(axes[0], order, young)
    panel_nano_sep(axes[1], young, adult)
    panel_swap(axes[2])

    fig.suptitle(
        "Choosing the reference channel: autofluorescence measures tissue, SEP was "
        "meant to measure the receptor and is mostly autofluorescence "
        "(run_sep_channel_check).\n"
        "Points off the identity line on the right are structures "
        "whose answer depends on that choice",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    save_figure(fig, OUT / "09_denominators.png", dpi=95, eps=False)
    plt.close(fig)


def write_index() -> None:
    """Write README.md beside the sheets: what to look for in each, what is wrong."""
    # the README as written; its table rows are longer than a line of code
    with open(OUT / "README.md", "w", encoding="utf-8") as fh:
        fh.write("""# processing diagnostics

One sheet per question, so every step of the v2 route can be checked by eye
rather than trusted. Regenerate with `run_diagnostics.py`.

| sheet | what to look for | what would be wrong |
|---|---|---|
| `01_tissue_<mouse>.png` | red contour hugs the tissue edge; ventricles and the space around fibre tracts excluded; blue atlas outline roughly matches | mask eating into cortex, or spilling into the black surround |
| `02_levels_<mouse>.png` | off-tissue and tissue distributions separate at the marked threshold; cortex centred on 1.0 after scaling | the two distributions overlapping at the threshold, or a cortex peak far from 1 |
| `03_coverage.png` | each brain covers a contiguous run of planes; the group overlaps in the middle | a brain with holes, or a cohort where only one brain covers a whole region |
| `04_warp_<mouse>.png` | before and after look like the same brain; region means sit on the identity line | points off the line, i.e. the transform moving signal between regions |
| `05_cohort_n.png` | the n map is flat across most of the brain | large areas resting on one or two brains |
| `06_scaling.png` | backgrounds similar across brains; cortex means spread widely; the two channels track each other | a background far from the others (a brain imaged differently) |
| `07_route_agreement.png` | points on the identity line | a systematic offset, i.e. the warp biasing the comparison |
| `08_mask_vs_p6bis.png` | red and blue contours on top of each other | the new mask cutting into tissue, or reaching into the surround, where the old one does not |
| `09_denominators.png` | the two reference channels behave alike across brains; structures sit on the identity line | a structure whose young-adult difference flips sign with the denominator -- that result belongs to the denominator, not to the biology |
| `sep_channel/<mouse>.png` | written by `run_add_sep_channel.m`: every slice recovered at r = 1.00000, and the registered DAPI of that run on top of the one run_register_to_atlas wrote | a slice below r = 0.99, or a residual shift above a fraction of a pixel: the SEP channel would not be sitting where NANO and AUTO sit |

The numbers behind these are in `../young_vs_adult/region_stats.csv` and
`../README.md`.
""")  # noqa: E501


def main(named_mice: list[str]) -> None:
    """Draw the sheets of `named_mice`, or of every brain and the cohort-level ones.

    With mouse names only those brains' sheets are refreshed; with none, every
    brain's sheets, the cohort-level sheets and the index are written. ISO and
    ANN, which the sheet functions read as module globals, are set here.
    """
    global ISO, ANN
    OUT.mkdir(parents=True, exist_ok=True)

    # the isocortex labels, and every atlas the brains are on
    ISO = list(isocortex_ids())
    ANN = {k: annotation_20(k) for k in {v[1] for v in MICE.values()}}

    # per-brain sheets; the warp sheet only for a young brain
    mice = named_mice or list(MICE)
    for mouse in mice:
        ann = ANN[MICE[mouse][1]]
        z = np.load(PER_MOUSE / (mouse + ".npz"))
        sheet_tissue(mouse, ann, z)
        sheet_levels(mouse, ann, z)
        if MICE[mouse][1] != "ccf":
            sheet_warp(mouse)
        print(f"{mouse:20s} sheets written", flush=True)

    # cohort-level sheets and the index, in a full run only
    if not named_mice:
        sheet_coverage()
        sheet_cohort_n()
        sheet_scaling()
        sheet_route_agreement()
        sheet_mask_vs_p6bis()
        sheet_denominators()
        write_index()
        print("cohort-level sheets and index written")
