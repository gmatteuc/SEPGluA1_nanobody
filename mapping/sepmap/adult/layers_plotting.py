"""The figures of the adult map by cortical depth (adult.layers).

    01  flatmaps by depth band: the mean of the ten adults, and the SD between them
    02  every isocortex area per band, along the Harris 2019 hierarchy, a dot per adult
    03  laminar profiles of the areas the argument is about, a line per adult

The flatmaps reuse the close-up's drawing (young_vs_adult.closeup.draw_flat: the
area borders, the names, the highlighted areas, black ground), so they read as the
adult column of the young against adult figure. Both rows use hot, cut at 0.82 of
its range, as every intensity map of the route: the means from adult_layers.mean_lim,
whose floor sits below the lowest band means so the darkest cortex stays dark red
rather than the black of the ground, and the SD from 0 to adult_layers.sd_max.

In the bar figure a bar is the mean over the adults with a value, its whisker the
SEM, and its grey the reliability of that mean, |t| = |mean| / SEM clamped at
adult_layers.t_max, darker for more reliable and never white, as the project's bars
are (sep_palette('bars')). The dots are the adults, in the nano channel's dot colour,
circles naive and triangles RWS.

Run by run_adult_layers.py.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from sepmap.config import SETTINGS
from sepmap.plotting import RED, hot_cut, save_figure, tidy
from sepmap.young_vs_adult import closeup

ADULT_LAYERS = SETTINGS["adult_layers"]

# the nano channel and its per-mouse dots, as sep_palette('nano') and ('nano_dot');
# plotting.py holds no channel colours yet
NANO = (0.95, 0.55, 0.10)
NANO_DOT = (0.65, 0.30, 0.00)

# the marker of each adult group
GROUP_MARKERS = {"naive": "o", "rws": "^"}

# the depth bands with the layers each pools, as the axis and titles name them
BAND_TITLES = {
    "supragranular": "supragranular  (L1 + L2/3)",
    "granular": "granular  (L4)",
    "infragranular": "infragranular  (L5 + L6a + L6b)",
}
BAND_ORDER = ["supragranular", "granular", "infragranular"]

# the five layers of a profile, top to bottom
PROFILE_LAYERS = ["L1", "L2/3", "L4", "L5", "L6"]

# the areas of the laminar profiles: the primary sensory areas, motor, retrosplenial
# and prefrontal cortex, then the three higher visual areas the grant is about
PROFILE_AREAS = [
    "SSp-bfd",
    "VISp",
    "AUDp",
    "MOp",
    "RSPd",
    "ACAd",
    "PL",
    "ORBl",
    "VISrl",
    "VISal",
    "VISli",
]

# the grant's areas, named in red
GRANT_AREAS = ("VISrl", "VISal", "VISli")


def bars_grey(t: np.ndarray, t_max: float) -> np.ndarray:
    """The grey of a bar of reliability `t`: 0.78 at 0, black at `t_max`, (n, 3) RGB.

    t is clamped to 0..t_max and NaN reads as 0, as sep_palette('bars') in MATLAB.
    """
    t = np.atleast_1d(np.asarray(t, dtype=float))
    share = np.clip(np.nan_to_num(t, nan=0.0) / t_max, 0.0, 1.0)
    level = 0.78 * (1.0 - share)
    return np.repeat(level[:, None], 3, axis=1)


def smoothing_text(sigma: float | list[float]) -> str:
    """The smoothing as the close-up's titles state it, in um along (AP, DV, ML)."""
    sig = np.atleast_1d(np.asarray(sigma, dtype=float))
    if not np.any(sig):
        return "no smoothing"
    sigmas_3d = sig if sig.size > 1 else np.repeat(sig, 3)
    return (
        "smoothed sigma "
        + " x ".join(f"{v * 20:.0f}" for v in sigmas_3d)
        + " um (AP x DV x ML)"
    )


def plot_band_flatmaps(
    mean: np.ndarray,
    sd: np.ndarray,
    n_adults: int,
    min_n: int,
    sigma: float | list[float],
    save: Path | None = None,
) -> plt.Figure:
    """Figure 01: the mean (top) and between-mouse SD (bottom) of each depth band.

    `mean` and `sd` are (band, row, col) in the order of closeup.BANDS; `min_n` is
    the adults a pixel needed.
    """
    left, right, label_xy = closeup.flatmap_borders()
    hot = hot_cut()
    fig, axes = plt.subplots(2, 3, figsize=(19, 9), facecolor="k")
    rows = (
        (mean, tuple(ADULT_LAYERS["mean_lim"]), "mean"),
        (sd, (0.0, ADULT_LAYERS["sd_max"]), "SD between adults"),
    )
    for ax_row, (maps, lim, what) in zip(axes, rows):
        for ax, img, (bname, _) in zip(ax_row, maps, closeup.BANDS):
            h = closeup.draw_flat(
                ax, img, hot, lim, f"{what}   {bname}", (left, right), label_xy
            )
            cb = fig.colorbar(h, ax=ax, fraction=0.03, pad=0.01)
            cb.ax.yaxis.set_tick_params(color="w", labelcolor="w")
    fig.suptitle(
        f"Adult isocortex by depth band, the {n_adults} adults: the mean of their band "
        f"maps and the SD between them  -  zref, {smoothing_text(sigma)}, pixels with "
        f"at least {min_n} adults",
        color="w",
        fontsize=12,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    if save is not None:
        save_figure(fig, save, dpi=150, facecolor="k")
    return fig


def area_positions(hier: pd.DataFrame) -> dict[str, float]:
    """x of each area: the scored ones by score, then a gap, then the unscored."""
    scored = hier["hierarchy_score"].notna()
    x = hier["order"] + (~scored).astype(int)
    return dict(zip(hier["acronym"], x.astype(float)))


def grey_key(ax: plt.Axes, t_max: float) -> None:
    """A small key to the grey of the bars, inside the panel at its upper right."""
    key = ax.inset_axes([0.80, 0.86, 0.16, 0.045])
    steps = np.linspace(0, t_max, 6)
    key.imshow(
        bars_grey(steps, t_max)[None, :, :], aspect="auto", extent=(0, t_max, 0, 1)
    )
    key.set_yticks([])
    key.set_xticks([0, t_max])
    key.set_xticklabels(["0", f"{t_max:g}"], fontsize=6.5)
    key.tick_params(length=2, pad=1)
    key.set_title("|mean| / SEM across the adults", fontsize=6.5, pad=2)


def draw_area_band(
    ax: plt.Axes,
    cells: pd.DataFrame,
    stats: pd.DataFrame,
    x_of: dict[str, float],
    t_max: float,
    rng: np.random.Generator,
) -> None:
    """One band's panel: a bar per area with its SEM, and a dot per adult."""
    # the bars: the mean, grey by its |t|, with the SEM
    s = stats[stats["mean"].notna()]
    xs = s["area"].map(x_of).to_numpy()
    ax.bar(
        xs,
        s["mean"],
        width=0.72,
        color=bars_grey(s["t"].abs().to_numpy(), t_max),
        edgecolor="none",
        zorder=1,
    )
    ax.errorbar(
        xs, s["mean"], yerr=s["sem"], fmt="none", ecolor="k", elinewidth=0.9, zorder=3
    )

    # the adults, jittered across the bar, by group
    kept = cells[~cells["excluded"]]
    for group, marker in GROUP_MARKERS.items():
        g = kept[kept["group"] == group]
        jitter = rng.uniform(-0.22, 0.22, len(g))
        ax.scatter(
            g["area"].map(x_of).to_numpy() + jitter,
            g["zref"],
            s=10,
            marker=marker,
            color=NANO_DOT,
            alpha=0.85,
            linewidths=0,
            zorder=4,
        )
    ax.axhline(0, color="k", lw=0.6, zorder=2)
    ax.grid(axis="y", lw=0.3, alpha=0.6)
    tidy(ax)


def band_title(band: str, depth_row: pd.Series, n_without: int) -> str:
    """A band panel's title: what it pools, its reliability and its hierarchy rho."""
    text = (
        f"{BAND_TITLES[band]}:   half against half rho {depth_row['half_rho']:.2f} over "
        f"the {depth_row['n_areas_complete']} areas every adult has, all ten "
        f"{depth_row['whole_rho']:.2f};   Spearman with the hierarchy "
        f"{depth_row['hierarchy_rho']:+.2f} (p = {depth_row['hierarchy_p']:.3f}, "
        f"{depth_row['n_areas_scored']} areas)"
    )
    if n_without:
        text += f";   {n_without} areas without layer 4 left blank"
    return text


def plot_area_bands(
    table: pd.DataFrame,
    summary: pd.DataFrame,
    depths: pd.DataFrame,
    hier: pd.DataFrame,
    n_adults: int,
    save: Path | None = None,
) -> plt.Figure:
    """Figure 02: every isocortex area per depth band, a dot per adult, mean and SEM.

    `table`, `summary` and `depths` are adult.layers' per-mouse table, area summary
    and depth summary; `hier` its hierarchy table.
    """
    x_of = area_positions(hier)
    n_scored = int(hier["hierarchy_score"].notna().sum())
    t_max = ADULT_LAYERS["t_max"]
    rng = np.random.default_rng(0)
    fig, axes = plt.subplots(3, 1, figsize=(17, 11.5), sharex=True)
    for ax, band in zip(axes, BAND_ORDER):
        cells = table[(table["depth_kind"] == "band") & (table["depth"] == band)]
        stats = summary[(summary["depth_kind"] == "band") & (summary["depth"] == band)]
        draw_area_band(ax, cells, stats, x_of, t_max, rng)
        row = depths[(depths["depth_kind"] == "band") & (depths["depth"] == band)].iloc[0]
        n_without = len(hier) - cells["area"].nunique()
        ax.set_title(band_title(band, row, n_without), fontsize=9, loc="left")
        ax.set_ylabel("zref")

        # the divide between the areas the hierarchy scores and those it does not
        ax.axvline(n_scored, color="k", lw=0.6, ls=":")
    grey_key(axes[0], t_max)

    # the areas named under the bottom panel, the grant's in red
    names = sorted(x_of, key=x_of.get)
    axes[-1].set_xticks([x_of[a] for a in names])
    axes[-1].set_xticklabels(names, rotation=90, fontsize=8)
    for label in axes[-1].get_xticklabels():
        if label.get_text() in GRANT_AREAS:
            label.set_color(RED)
    axes[-1].set_xlim(-0.8, max(x_of.values()) + 0.8)
    axes[-1].set_xlabel(
        f"isocortex areas, from low to high in the hierarchy of Harris et al. 2019 "
        f"({n_scored} areas, Supplementary Table 9); right of the dotted line, the "
        f"{len(hier) - n_scored} it does not score"
    )

    # the legend of the dots and the bars, once, on the top panel
    handles = [
        Line2D([], [], marker=m, ls="none", color=NANO_DOT, label=f"adult {g}")
        for g, m in GROUP_MARKERS.items()
    ]
    handles.append(Line2D([], [], color="k", lw=0.9, label="mean and SEM"))
    axes[0].legend(handles=handles, loc="upper left", fontsize=8, frameon=False, ncol=3)
    fig.suptitle(
        f"The adult isocortex by depth band, area by area: one dot per adult (n = "
        f"{n_adults}), the bar their mean, grey by how reliable it is\n"
        "zref: 0 is the brain's median structure, 1 its p90 - p10 spread across "
        "structures",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    if save is not None:
        save_figure(fig, save, dpi=150)
    return fig


def draw_profile(
    ax: plt.Axes, cells: pd.DataFrame, stats: pd.DataFrame, title: str, colour: str
) -> None:
    """One area's laminar profile: a line per adult, the mean with its SEM on top."""
    # the layers the area has, top to bottom; the others are skipped by the lines
    present = [lay for lay in PROFILE_LAYERS if lay in set(stats["depth"])]
    y_of = {lay: PROFILE_LAYERS.index(lay) for lay in PROFILE_LAYERS}
    for (mouse, group), part in cells.groupby(["mouse", "group"], sort=False):
        part = part.set_index("depth").reindex(present)
        ax.plot(
            part["zref"],
            [y_of[lay] for lay in present],
            color=NANO_DOT,
            lw=0.8,
            alpha=0.55,
            ls="-" if group == "naive" else "--",
            marker=GROUP_MARKERS[group],
            ms=2.5,
        )
    s = stats.set_index("depth").reindex(present)
    ys = [y_of[lay] for lay in present]
    ax.errorbar(
        s["mean"], ys, xerr=s["sem"], color=NANO, lw=2.2, elinewidth=1.2, zorder=5
    )
    if "L4" not in present:
        ax.text(
            0.03,
            y_of["L4"],
            "no layer 4",
            ha="left",
            va="center",
            fontsize=7,
            color="0.5",
            transform=ax.get_yaxis_transform(),
        )
    ax.axvline(0, color="k", lw=0.6)
    ax.set_yticks(range(len(PROFILE_LAYERS)))
    ax.set_yticklabels(PROFILE_LAYERS)
    ax.set_ylim(len(PROFILE_LAYERS) - 0.5, -0.5)
    ax.set_title(title, fontsize=9, color=colour)
    ax.grid(axis="x", lw=0.3, alpha=0.6)
    tidy(ax)


def profile_title(area: str, hier: pd.DataFrame) -> str:
    """An area's panel title: its module and its rank in the hierarchy."""
    row = hier.set_index("acronym").loc[area]
    n_scored = int(hier["hierarchy_score"].notna().sum())
    if np.isnan(row["hierarchy_rank"]):
        rank = "not scored"
    else:
        rank = f"{int(row['hierarchy_rank'])} of {n_scored}"
    return f"{area}   {row['module'].lower()}, hierarchy {rank}"


def plot_laminar_profiles(
    table: pd.DataFrame,
    summary: pd.DataFrame,
    hier: pd.DataFrame,
    n_adults: int,
    save: Path | None = None,
) -> plt.Figure:
    """Figure 03: the five-layer profile of each of PROFILE_AREAS, a line per adult."""
    cells = table[table["depth_kind"] == "layer"]
    stats = summary[summary["depth_kind"] == "layer"]
    shown = cells[cells["area"].isin(PROFILE_AREAS) & ~cells["excluded"]]

    # one x range for every panel, so the areas compare by eye
    lo, hi = np.nanmin(shown["zref"]), np.nanmax(shown["zref"])
    pad = 0.05 * (hi - lo)
    fig, axes = plt.subplots(3, 4, figsize=(14, 10), sharex=True)
    for ax, area in zip(axes.flat, PROFILE_AREAS):
        colour = RED if area in GRANT_AREAS else "k"
        draw_profile(
            ax,
            shown[shown["area"] == area],
            stats[stats["area"] == area],
            profile_title(area, hier),
            colour,
        )
        ax.set_xlim(lo - pad, hi + pad)

    # the last cell holds the legend
    key = axes.flat[-1]
    key.axis("off")
    handles = [
        Line2D(
            [],
            [],
            color=NANO_DOT,
            lw=0.8,
            ls="-",
            marker="o",
            ms=2.5,
            label="adult naive",
        ),
        Line2D(
            [], [], color=NANO_DOT, lw=0.8, ls="--", marker="^", ms=2.5, label="adult rws"
        ),
        Line2D([], [], color=NANO, lw=2.2, label="mean and SEM"),
    ]
    key.legend(handles=handles, loc="center", fontsize=9, frameon=False)

    # the x axis named under the lowest panel of each column, which the legend may be
    n_cols = axes.shape[1]
    for col in range(n_cols):
        lowest = max(i for i in range(len(PROFILE_AREAS)) if i % n_cols == col)
        axes.flat[lowest].tick_params(labelbottom=True)
        axes.flat[lowest].set_xlabel("zref")
    fig.suptitle(
        f"Laminar profiles of the adult map, one line per adult (n = {n_adults}); "
        "L6 pools 6a and 6b\nzref: 0 is the brain's median structure, 1 its p90 - p10 "
        "spread across structures",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    if save is not None:
        save_figure(fig, save, dpi=150)
    return fig
