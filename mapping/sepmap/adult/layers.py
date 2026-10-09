"""The adult map by cortical depth: every isocortex area, depth by depth, mouse by mouse.

The young against adult close-ups unroll the isocortex into flatmaps by depth band;
this module gives the adult map the same treatment on its own, and beside it the
numbers per mouse, so the variability from one adult to the next can be seen area
by area and layer by layer.

The reading is zref as the region tables compute it (young_vs_adult.region_plot):
per brain, the voxel-weighted mean nano (sig) of a set of labels, then

    zref = ( log2(mean / isocortex mean) - median ) / (p90 - p10)

with the median and the spread of that brain's structures taken over the structures
every brain of the region tables shares (region_plot.range_match, over the young
brains and the adults alike). The same function on the same brains, so a cell taken
over an area's whole depth equals that area's row of region_means_per_mouse.csv,
which the run checks. Zero is the brain's median structure, one unit its own
spread; no zref value is a fold change.

The cells, per adult, area and depth:

    areas    the 43 isocortex areas of the structure term set, as in the region
             tables
    bands    supragranular (layers 1, 2/3), granular (4) and infragranular (5, 6a,
             6b), the bands of young_vs_adult.region_groups and of the close-up
    layers   L1, L2/3, L4, L5 and L6 (6a and 6b pooled), for the laminar profiles
    whole    the area through its whole depth, for the check
    contrast supragranular minus infragranular, mouse by mouse: a difference
             within one brain, so its level owes nothing to that brain's median

A cell is the voxel-weighted mean over the area's labels of that depth in that
brain (region_groups.group_cells), kept from region_tables.min_vox20 tissue voxels
(250, 2 nl), as a structure is in the region tables. A smaller cell, mostly an area
at the edge of what the sections reached, is missing, and the table says why. The
CCF draws no voxel of RSPd layer 4 though the ontology names it, so a depth with no
voxel in the atlas is left out of every table. Across the adults, an area and depth
gets a mean, SD, SEM and t = mean / SEM from adult_layers.min_mice adults.

The order of the areas is the cortical hierarchy of Harris et al. 2019 (Nature
575:195), the final hierarchy of their Fig. 6d from corticocortical,
thalamocortical and corticothalamic connections with Cre-line confidence: column
"Cre_conf CC+TC+CT" of the tab "all area hr scores_Cre + WT" of Supplementary Table
9 (41586_2019_1716_MOESM10_ESM.xlsx, MD5 55b69aeba95468559ff025e3c08b24b6, the MD5
that PubMed Central lists for the same table), fetched on 9 October 2026. All 37
scores equal those of the authors' code, github.com/AllenInstitute/MouseBrainHierarchy
at commit 8e4e0dd (Results/hierarchy_summary_CreConf.xlsx, "CC+TC+CT iterated", to
1e-9; Output/TCCT_CCconf_iter.xls to 1e-16), and as the paper says, VISp is lowest
and ORBvl highest. The modules are the authors' too (Input/CC_TC_CT_clusters.xlsx),
and they name all 43 areas. Six areas carry no score, AUDv, SSp-un, ECT, GU, PERI
and VISC; they follow the 37, in the order of Harris's hierarchy of modules (Fig.
6f: auditory, visual, somatomotor, medial, lateral, prefrontal). The table, with
the CCF's spelling of AId, AIp and AIv, is adult_layers.hierarchy.

Two numbers per depth say how far the profile across areas can be trusted and how
it lies along the hierarchy:

    reliability   the ten adults cut into two fives every possible way (126 cuts),
                  the two half-cohort profiles over the areas every adult has
                  compared by Spearman, the mean of that, and the Spearman-Brown
                  value for all ten, 2 r / (1 + r), as adult.beyond_density takes
                  the ceiling of the whole map
    hierarchy     Spearman of the area means with the hierarchy scores. Descriptive:
                  neighbouring areas share signal and a hierarchy, so its p treats
                  as independent what is not

The flatmaps are drawn from the voxel maps, as the close-up's: each adult's zref
volume (volumes.cohort.mouse_modes, whose median and spread come from the brain's
own structures, which is why the route quotes numbers from its tables and reads
pattern from its maps), hemispheres folded, smoothed inside its tissue as the
close-up smooths a cohort (closeup.smooth_within, closeup.smooth), projected along
the streamlines into depth bins (closeup.project_slab) and averaged over each band
(closeup.band_average). Per pixel, the mean and the SD of the ten band maps, from
young_vs_adult.min_n_adult adults, the threshold of the close-up's adult panel.
Taking the SD over each adult's own band map makes it the variability between
mice of what the mean shows, not that of single voxels. The projection takes about
a minute per brain, so the ten band maps are cached.

Writes, in adult_v2/layers/ under the data root:

    area_layers_per_mouse.csv  one row per adult, area and depth, kept or missing
                               with the reason
    area_layers_summary.csv    one row per area and depth: n, mean, SD, SEM, t, the
                               adults missing
    depth_summary.csv          per depth, the reliability and the hierarchy rho
    flatmap_bands.npz          the cached band maps of every adult

and the figures of adult.layers_plotting.

Run by run_adult_layers.py.
"""

import itertools
import re
import warnings
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from sepmap.config import DATA, SETTINGS, code_root
from sepmap.hemispheres import fold
from sepmap.volumes.cohort import NAIVE, RWS, mouse_modes
from sepmap.volumes.per_mouse import MICE, annotation_20
from sepmap.volumes.per_mouse import OUT as PER_MOUSE
from sepmap.young_vs_adult import closeup
from sepmap.young_vs_adult.region_groups import (
    LAYERS,
    group_cells,
    layer_of,
    references,
    structure_sig_means,
)
from sepmap.young_vs_adult.region_plot import GROUPS, REGION_MEANS, label_sums, value

ADULT_LAYERS = SETTINGS["adult_layers"]
REGION_TABLES = SETTINGS["region_tables"]

OUT = DATA / "adult_v2" / "layers"
PER_MOUSE_CCF = DATA / "comparisons_v2" / "per_mouse_ccf"

# the ten adults, naive and rws pooled, and every brain of the region tables, whose
# shared structures set each brain's zref median and spread
ADULTS = NAIVE + RWS
ROUTE_MICE = [m for g in GROUPS.values() for m in g]

# the three depth bands, in the order of young_vs_adult.region_groups
BANDS = [name for name, _ in LAYERS]

# the five layers of the laminar profiles: layer 6 pools 6a and 6b
LAYERS_5 = [
    ("L1", ("1",)),
    ("L2/3", ("2/3", "2", "3")),
    ("L4", ("4",)),
    ("L5", ("5",)),
    ("L6", ("6", "6a", "6b")),
]

# the laminar contrast, taken mouse by mouse
CONTRAST = "supragranular - infragranular"

# the depths in the order of the tables: the whole depth, the bands, the layers, the
# contrast
DEPTHS = (
    [("whole", "whole")]
    + [("band", name) for name in BANDS]
    + [("layer", name) for name, _ in LAYERS_5]
    + [("contrast", CONTRAST)]
)

# the columns of the two tables
PER_MOUSE_COLUMNS = [
    "reading",
    "depth_kind",
    "depth",
    "area",
    "module",
    "hierarchy_score",
    "hierarchy_rank",
    "mouse",
    "group",
    "n_vox20",
    "atlas_vox20",
    "coverage",
    "zref",
    "excluded",
    "exclude_reason",
]
SUMMARY_COLUMNS = [
    "depth_kind",
    "depth",
    "area",
    "module",
    "hierarchy_score",
    "hierarchy_rank",
    "n_mice",
    "mean",
    "sd",
    "sem",
    "t",
    "mice_missing",
]


def load_hierarchy(path: str | Path | None = None) -> pd.DataFrame:
    """The hierarchy table: acronym, module, hierarchy_score, rank and plotting order.

    `path` defaults to adult_layers.hierarchy under the code root. hierarchy_rank
    runs from 1 (lowest) over the scored areas and is NaN for the others; `order`
    is the position along a figure's axis, the scored areas by score, then the
    others in the file's order. Stops on a repeated acronym.
    """
    if path is None:
        path = Path(code_root()) / ADULT_LAYERS["hierarchy"]
    hier = pd.read_csv(path)
    if hier["acronym"].duplicated().any():
        repeated = hier.loc[hier["acronym"].duplicated(), "acronym"].tolist()
        raise ValueError(f"{path}: acronyms listed twice: {repeated}")

    # rank over the scored areas, then the unscored after them, as listed
    hier["hierarchy_rank"] = hier["hierarchy_score"].rank(method="first")
    scored = hier["hierarchy_score"].notna()
    order = pd.Series(np.nan, index=hier.index)
    order[scored] = hier.loc[scored, "hierarchy_rank"] - 1
    order[~scored] = scored.sum() + np.arange((~scored).sum())
    hier["order"] = order.astype(int)
    return hier.sort_values("order").reset_index(drop=True)


def layer_token(substructure_name: str) -> str | None:
    """The layer of an isocortex label ('2/3', '6a'), or None without one.

    As region_groups.layer_of reads it, but for two labels that give the layer
    without the word: "Anterior cingulate area, ventral part, 6a" and "..., 6b".
    layer_of returns None for them, so the layer is then the last comma-separated
    token, when it is a layer.
    """
    token = layer_of(substructure_name)
    if token is not None:
        return token
    last = substructure_name.rsplit(",", 1)[-1].strip()
    if re.fullmatch(r"[0-9]+(?:/[0-9]+)?[ab]?", last):
        return last
    return None


def depth_groups(
    stru: dict[int, str], divi: dict[int, str], sub: dict[int, str]
) -> dict[tuple[str, str, str], set[int]]:
    """The labels of every isocortex area at every depth, keyed (kind, depth, area).

    `stru`, `divi` and `sub` give the structure and division acronyms and the
    substructure name of each parcellation index (region_groups'
    load_parcellation_terms). Kinds: "band" (BANDS), "layer" (LAYERS_5) and
    "whole". A depth an area has no label of is not a key.
    """
    groups = {}
    iso = [i for i in stru if divi.get(i) == "Isocortex"]
    token = {i: layer_token(sub.get(i, "")) for i in iso}
    for area in sorted({stru[i] for i in iso}):
        ids = [i for i in iso if stru[i] == area]
        groups[("whole", "whole", area)] = set(ids)
        for name, tokens in LAYERS:
            groups[("band", name, area)] = {i for i in ids if token[i] in tokens}
        for name, tokens in LAYERS_5:
            groups[("layer", name, area)] = {i for i in ids if token[i] in tokens}
    return {k: v for k, v in groups.items() if v}


def check_hierarchy(hier: pd.DataFrame, areas: set[str]) -> None:
    """Stop unless the hierarchy table names exactly the isocortex areas of the atlas."""
    listed = set(hier["acronym"])
    if listed != areas:
        raise ValueError(
            f"the hierarchy table ({ADULT_LAYERS['hierarchy']}) does not match the "
            f"atlas's isocortex areas: missing {sorted(areas - listed)}, not in the "
            f"atlas {sorted(listed - areas)}"
        )


def atlas_counts(groups: dict[tuple[str, str, str], set[int]]) -> dict[tuple, int]:
    """Voxels of each group in the adults' atlas, the cropped CCF at 20 um."""
    ann = annotation_20("ccf")
    n = np.bincount(ann.ravel(), minlength=int(ann.max()) + 1)
    return {k: int(n[[i for i in ids if i < len(n)]].sum()) for k, ids in groups.items()}


def measure(
    groups: dict[tuple[str, str, str], set[int]],
    stru: dict[int, str],
    divi: dict[int, str],
) -> tuple[dict, dict, dict, dict]:
    """Per brain: the cells of the adults, their voxel counts, the references, the
    structure means.

    Every brain of the region tables (ROUTE_MICE) gives its isocortex and subcortex
    means and its structure means, which set the zref of every brain; the adults
    also give a cell per group (region_groups.group_cells: voxels, mean sig, mean
    ratio, mean sepratio, None under region_tables.min_vox20 voxels) and the tissue
    voxels of each group, kept for the cells that fall under it. Returns (cells,
    n_vox, refs, struct_mean), each keyed by mouse.
    """
    anns, cells, n_vox, refs, struct_mean = {}, {}, {}, {}, {}
    for mouse in ROUTE_MICE:
        # the brain's own atlas, loaded once for all the brains on it
        atlas_key = MICE[mouse][1]
        if atlas_key not in anns:
            anns[atlas_key] = annotation_20(atlas_key)

        # sums per parcellation index, the two references and the structure means
        sums = label_sums(mouse, anns[atlas_key])
        n, s_sig = sums["n"], sums["sig"]
        refs[mouse] = references(stru, divi, n, s_sig, len(n))
        struct_mean[mouse] = structure_sig_means(stru, n, s_sig)

        # the groups, for the adults only
        if mouse in ADULTS:
            cells[mouse] = group_cells(groups, sums)
            n_vox[mouse] = {
                k: int(n[[i for i in ids if i < len(n)]].sum())
                for k, ids in groups.items()
            }
            n_kept = sum(c is not None for c in cells[mouse].values())
            print(f"  {mouse:20s} {n_kept}/{len(groups)} cells kept", flush=True)
    return cells, n_vox, refs, struct_mean


def exclusion(n: int, zref: float | None) -> str:
    """Why a cell has no value, or '' when it has one."""
    if zref is not None:
        return ""
    if n == 0:
        return "no tissue in the sections"
    if n < REGION_TABLES["min_vox20"]:
        return f"{n} tissue voxels, under region_tables.min_vox20"
    return "mean nano at or below background"


def _mouse_rows(
    key: tuple[str, str, str],
    cells: dict,
    n_vox: dict,
    norm: dict[str, tuple[float, float]],
    refs: dict[str, dict[str, float]],
    by_area: pd.DataFrame,
    atlas_n: dict[tuple, int],
) -> list[dict]:
    """The rows of one group, one per adult."""
    kind, depth, area = key
    rows = []
    for mouse in ADULTS:
        z = value("zref", mouse, key, cells, norm, refs)
        rows.append(
            dict(
                reading="zref",
                depth_kind=kind,
                depth=depth,
                area=area,
                module=by_area.at[area, "module"],
                hierarchy_score=by_area.at[area, "hierarchy_score"],
                hierarchy_rank=by_area.at[area, "hierarchy_rank"],
                mouse=mouse,
                group=MICE[mouse][0],
                n_vox20=n_vox[mouse][key],
                atlas_vox20=atlas_n[key],
                coverage=n_vox[mouse][key] / atlas_n[key],
                zref=np.nan if z is None else z,
                excluded=z is None,
                exclude_reason=exclusion(n_vox[mouse][key], z),
            )
        )
    return rows


def per_mouse_table(
    groups: dict[tuple[str, str, str], set[int]],
    cells: dict,
    n_vox: dict,
    norm: dict[str, tuple[float, float]],
    refs: dict[str, dict[str, float]],
    hier: pd.DataFrame,
    atlas_n: dict[tuple, int],
) -> pd.DataFrame:
    """One row per adult and group the atlas draws, with its zref or why it has none.

    `norm` is region_plot.range_match's median and spread per brain. A group with
    no voxel in the atlas (RSPd layer 4) gets no row. Rows run by depth (DEPTHS),
    then along the hierarchy, then by adult.
    """
    by_area = hier.set_index("acronym")
    rows = []
    for kind, depth in DEPTHS:
        for area in hier["acronym"]:
            key = (kind, depth, area)
            if key in groups and atlas_n[key] > 0:
                rows.extend(_mouse_rows(key, cells, n_vox, norm, refs, by_area, atlas_n))
    return pd.DataFrame(rows, columns=PER_MOUSE_COLUMNS)


def with_contrast(table: pd.DataFrame) -> pd.DataFrame:
    """The per-mouse table with the supragranular minus infragranular rows added.

    Per adult and area, a difference of two cells of the same brain; missing when
    either is, with the reason of the first that is. Its voxel counts are the
    smaller of the two, its coverage too.
    """
    bands = table[table["depth_kind"] == "band"].set_index(["area", "mouse", "depth"])
    rows = []
    for (area, mouse), part in bands.groupby(level=["area", "mouse"], sort=False):
        part = part.droplevel(["area", "mouse"])
        if "supragranular" not in part.index or "infragranular" not in part.index:
            continue
        supra, infra = part.loc["supragranular"], part.loc["infragranular"]
        excluded = bool(supra["excluded"] or infra["excluded"])
        reason = ""
        if supra["excluded"]:
            reason = "supragranular: " + supra["exclude_reason"]
        elif infra["excluded"]:
            reason = "infragranular: " + infra["exclude_reason"]
        row = supra.to_dict()
        row.update(
            depth_kind="contrast",
            depth=CONTRAST,
            area=area,
            mouse=mouse,
            n_vox20=min(supra["n_vox20"], infra["n_vox20"]),
            atlas_vox20=min(supra["atlas_vox20"], infra["atlas_vox20"]),
            coverage=min(supra["coverage"], infra["coverage"]),
            zref=np.nan if excluded else supra["zref"] - infra["zref"],
            excluded=excluded,
            exclude_reason=reason,
        )
        rows.append(row)
    contrast = pd.DataFrame(rows, columns=PER_MOUSE_COLUMNS)
    return pd.concat([table, contrast], ignore_index=True)


def summary_table(table: pd.DataFrame, min_mice: int | None = None) -> pd.DataFrame:
    """One row per area and depth: the adults with a value, their mean, SD, SEM and t.

    t is mean / SEM, signed. The statistics are NaN with fewer than `min_mice` adults
    (adult_layers.min_mice); mice_missing names the adults without a value,
    semicolon separated.
    """
    if min_mice is None:
        min_mice = ADULT_LAYERS["min_mice"]
    rows = []
    keys = ["depth_kind", "depth", "area"]
    for (kind, depth, area), part in table.groupby(keys, sort=False):
        kept = part.loc[~part["excluded"], "zref"].to_numpy(dtype=float)
        n = len(kept)
        mean = sd = sem = t = np.nan
        if n >= min_mice:
            mean = float(kept.mean())
            sd = float(kept.std(ddof=1))
            sem = sd / np.sqrt(n)
            t = mean / sem if sem > 0 else np.nan
        first = part.iloc[0]
        rows.append(
            dict(
                depth_kind=kind,
                depth=depth,
                area=area,
                module=first["module"],
                hierarchy_score=first["hierarchy_score"],
                hierarchy_rank=first["hierarchy_rank"],
                n_mice=n,
                mean=mean,
                sd=sd,
                sem=sem,
                t=t,
                mice_missing=";".join(part.loc[part["excluded"], "mouse"]),
            )
        )
    return pd.DataFrame(rows, columns=SUMMARY_COLUMNS)


def split_half(values: np.ndarray) -> tuple[float, float, int]:
    """Agreement of two half-cohort profiles over every split, and its Spearman-Brown.

    `values` is (mice, areas), complete. Every way of cutting the mice into two
    halves is taken once (the halves holding mouse 0), each half's mean profile is
    compared with the other's by Spearman, and the mean over the cuts r gives the
    whole cohort's 2 r / (1 + r). Returns (r, whole, number of cuts).
    """
    n_mice = values.shape[0]
    rhos = []
    for first in itertools.combinations(range(n_mice), n_mice // 2):
        if 0 not in first:
            continue
        second = [i for i in range(n_mice) if i not in first]
        a = values[list(first)].mean(axis=0)
        b = values[second].mean(axis=0)
        rhos.append(spearmanr(a, b).statistic)
    r = float(np.mean(rhos))
    whole = 2 * r / (1 + r) if r > -1 else float("nan")
    return r, whole, len(rhos)


def depth_summary(table: pd.DataFrame, summary: pd.DataFrame) -> pd.DataFrame:
    """Per depth: the reliability of the profile across areas, and its hierarchy rho.

    The reliability is taken over the areas every adult has at that depth
    (split_half); the hierarchy rho over the scored areas with a group mean.
    """
    rows = []
    for (kind, depth), part in table.groupby(["depth_kind", "depth"], sort=False):
        # the areas every adult has, as a (mice, areas) array
        wide = part.pivot(index="mouse", columns="area", values="zref")
        complete = wide.dropna(axis=1)
        r = whole = np.nan
        n_cuts = 0
        if complete.shape[1] >= 3:
            r, whole, n_cuts = split_half(complete.to_numpy())

        # the area means against the hierarchy scores
        s = summary[(summary["depth_kind"] == kind) & (summary["depth"] == depth)]
        s = s[s["mean"].notna() & s["hierarchy_score"].notna()]
        rho = p = np.nan
        if len(s) >= 3:
            res = spearmanr(s["hierarchy_score"], s["mean"])
            rho, p = float(res.statistic), float(res.pvalue)
        rows.append(
            dict(
                depth_kind=kind,
                depth=depth,
                n_areas_complete=complete.shape[1],
                n_cuts=n_cuts,
                half_rho=r,
                whole_rho=whole,
                n_areas_scored=len(s),
                hierarchy_rho=rho,
                hierarchy_p=p,
            )
        )
    return pd.DataFrame(rows)


def check_against_region_table(table: pd.DataFrame, path: Path = REGION_MEANS) -> int:
    """Stop unless every whole-depth cell equals its row of the region table.

    The region table holds zref to four decimals, so a cell may differ by half the
    last of them; a larger gap, or an area kept here and missing there, means the
    two were made from different per-mouse files or brains. Returns the cells
    compared.
    """
    region = pd.read_csv(path)
    region = region[(region["reading"] == "zref") & (region["division"] == "Isocortex")]
    region = region.set_index(["acronym", "mouse"])["log2_value"]
    whole = table[(table["depth_kind"] == "whole") & ~table["excluded"]]
    whole = whole.set_index(["area", "mouse"])["zref"]
    missing = whole.index.difference(region.index)
    if len(missing):
        raise ValueError(
            f"{len(missing)} whole-area cells have no row in {path}, e.g. "
            f"{list(missing[:3])}; run run_region_plot.py again on these per-mouse files"
        )
    gap = (whole - region.loc[whole.index]).abs().max()
    if gap > 0.5e-4 + 1e-9:
        raise ValueError(
            f"whole-area zref differs from {path} by up to {gap:.5f}; "
            "run run_region_plot.py again on these per-mouse files"
        )
    return len(whole)


def source_dates(mice: Sequence[str]) -> np.ndarray:
    """The modification times of each adult's per-mouse files, CCF and own atlas."""
    return np.array(
        [
            [
                (PER_MOUSE_CCF / f"{m}.npz").stat().st_mtime,
                (PER_MOUSE / f"{m}.npz").stat().st_mtime,
            ]
            for m in mice
        ]
    )


def mouse_band_maps(mouse: str, p3, sigma: float | list[float]) -> np.ndarray:
    """One adult's zref averaged over each depth band of the close-up, as flatmaps.

    The voxel zref of the maps (volumes.cohort.mouse_modes), hemispheres folded and
    smoothed inside the tissue as the close-up prepares a cohort, projected into
    depth bins. Returns (band, row, col), NaN where a band holds none of the brain's
    tissue.
    """
    volumes, _ = mouse_modes(mouse)

    # a voxel with tissue on neither side stays NaN, which nanmean warns about
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        half = fold(volumes["zref"])
    del volumes
    mask = np.isfinite(half).astype(np.float32)
    v = np.where(mask > 0, half, 0).astype(np.float32)
    if np.any(sigma):
        v = closeup.smooth_within(v, mask, sigma)

    # the slab of depth bins, then the mean over each band's bins
    sv, sm = closeup.project_slab(closeup.to_10um(v), closeup.to_10um(mask), p3)
    edges = closeup.band_edges(p3, sv.shape[2])
    maps = []
    for _, keys in closeup.BANDS:
        lo, hi = edges[keys[0]][0], edges[keys[-1]][1]
        maps.append(closeup.band_average(sv, sm, lo, hi))
    return np.stack(maps).astype(np.float32)


def band_maps(sigma: float | list[float], recompute: bool = False) -> np.ndarray:
    """The band flatmaps of every adult, (mouse, band, row, col), cached.

    The cache (adult_layers.flatmap_cache under OUT) is used when it holds the same
    adults, the same smoothing and per-mouse files of the same dates; `recompute`
    projects again whatever it holds.
    """
    cache = OUT / ADULT_LAYERS["flatmap_cache"]
    sig = np.atleast_1d(np.asarray(sigma, dtype=float))
    dates = source_dates(ADULTS)
    if cache.exists() and not recompute:
        z = np.load(cache)
        same = (
            list(z["mice"]) == ADULTS
            and np.array_equal(z["sigma"], sig)
            and np.array_equal(z["dates"], dates)
        )
        if same:
            print(
                f"  band maps of {len(ADULTS)} adults read from {cache.name}", flush=True
            )
            return z["maps"]
        print(f"  {cache.name} was made from other inputs; projecting again", flush=True)

    # the projector into depth bins, once for every brain
    _, p3 = closeup.flatmap_projectors()
    maps = []
    for mouse in ADULTS:
        maps.append(mouse_band_maps(mouse, p3, sigma))
        print(f"  projected {mouse}", flush=True)
    maps = np.stack(maps)
    np.savez_compressed(cache, maps=maps, mice=np.array(ADULTS), sigma=sig, dates=dates)
    return maps


def band_stats(maps: np.ndarray, min_n: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per pixel and band, the mean, SD and count across the mice of (mouse, band, ...).

    Mean and SD are NaN under `min_n` mice with a value, the SD (with n - 1) also
    under two.
    """
    n = np.isfinite(maps).sum(axis=0)

    # a pixel no mouse reaches is all NaN, which nanmean and nanstd warn about
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        mean = np.nanmean(maps, axis=0)
        sd = np.nanstd(maps, axis=0, ddof=1)
    keep = n >= max(min_n, 2)
    return np.where(keep, mean, np.nan), np.where(keep, sd, np.nan), n
