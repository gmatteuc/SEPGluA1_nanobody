"""Cohort volumes in the adult CCF: mean, SD and n per voxel, for every reading.

Everything is averaged on one grid, 660 x 400 x 570 at 20 um, where
volumes.to_ccf has carried each brain: the young through the DeMBA -> CCF
transform of their own age, the adults by placement alone. A pooled young group
therefore mixes ages without one age being carried by another's deformation
field.

Five readings of the same background-subtracted signal, the five the region
tables carry, so a number in a table and a colour in a map mean the same thing.
Two divide by a channel, voxel by voxel; three by a number measured on the
brain itself:

    ratio     sig / auto per voxel (auto smoothed by one 20 um voxel, so a dark
              voxel cannot blow it up), nano per unit autofluorescence. Not an
              absolute measurement: in the isocortex the young sit 2.0 log2 below
              the adults in nano and 1.0 log2 below them in auto, and
              autofluorescence rises with age (lipofuscin, tissue density), so
              this reading understates a real pup deficit and would overstate a
              pup excess. A bound, not a value.
    sepratio  sig / SEP per voxel, meant as surface receptor per unit receptor
              expressed. The green channel in this fixed, mounted tissue is
              mostly autofluorescence (rho 0.79 +- 0.04 against the auto channel
              in all ten adults; dynamic range 0.95 log2 against auto's 1.07 and
              nano's 1.93; adult.sep_channel_check), so this reading tracks ratio
              at rho 0.89 to 0.97 within every mouse. It is not a surface fraction.
    cref      sig / that mouse's isocortex mean of sig, measured before any warp:
              a pure scale, so region ratios within a mouse survive exactly; blind
              to a change that moves the whole cortex.
    subref    sig / that mouse's subcortex mean: TH, HY, PAL, MB, P and MY, the
              divisions not in NOT_SUBCORTEX (which holds the isocortex, OLF, HPF,
              CTXsp, STR, CB, the fiber tracts, the ventricles and the atlas's two
              catch-all labels). A reference that is neither the cortex nor HPF and
              STR, which would dominate the scale, and holds no white matter,
              ventricle or voxel of no structure.
    zref      log2(sig / cortex mean), minus that brain's median over structures,
              divided by its own p90-p10 spread, so level and dynamic range are
              the same in every brain. Median and spread are taken over the
              brain's structures of at least region_tables.min_vox20 voxels
              (mouse_scalars). The pup brain is flatter: its spread is about half
              the adult one, 0.9 to 1.0 log2 against 1.8 to 1.9, depending on
              which brains and structures enter. The only reading not linear in
              the signal: its maps show a position within a range, and the
              compression it removes may itself be the finding, so it is read
              beside cref.

Per voxel a cohort gets the mean over the mice with a value there (tissue, and for
ratio and sepratio a positive reference), the SD and that count. No voxel needs
every mouse; the n map says what each mean rests on, and the reader thresholds
it. Cohorts: young (P16, P20 and P22 pooled, the main
comparison), young_P20 (the P20 brains alone, the sensitivity check), young_P16,
young_P22, naive, rws, and adult (naive and rws).

The same statistics are taken a second time with each brain's two hemispheres
averaged first (hemispheres.fold), on the left half of the grid.
Each brain then gives one value per voxel, so the count is the brains with a value
on either side and a t over them has the n it claims; folding the cohort's own
left and right maps instead leaves no such n where the two sides were cut in
different brains.

Writes comparisons_v2/ccf/<cohort>/<reading>_{mean,sd,n}.npy,
<reading>_folded_{mean,sd,n}.npy (660 x 400 x 285) and mice.txt under the data
root.

Run by run_cohort.py; its cohorts and readings are imported across the package.
"""

import os
import warnings
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter

from sepmap.config import SETTINGS
from sepmap.hemispheres import fold
from sepmap.volumes.per_mouse import DATA, MICE, annotation_20, structure_terms
from sepmap.volumes.per_mouse import OUT as PER_MOUSE

READINGS = SETTINGS["readings"]
REGION_TABLES = SETTINGS["region_tables"]

V2 = DATA / "comparisons_v2"
PER_MOUSE_CCF = V2 / "per_mouse_ccf"
OUT_ROOT = V2 / "ccf"

# every reading the route can compute, in the order of the tables
ALL_MODES = ("ratio", "sepratio", "cref", "subref", "zref")

# the readings that are a signed position rather than an intensity: young and adult
# are compared by difference, not log2 ratio, on a diverging scale centred on zero.
# Defined once, in settings.toml, because one script getting it wrong would look like
# a result, not a bug; the chain imports it from here, and young_vs_adult.closeup
# reads the same key.
SIGNED_READINGS = tuple(READINGS["signed"])

# the readings in force, readings.in_force of settings.toml, or for one run
# V2_READINGS, a comma-separated subset ($env:V2_READINGS = 'ratio,cref' in
# PowerShell). The others are dropped from every figure, table and video without
# touching the data; young_vs_adult.region_plot and region_groups follow it too.
_want = os.environ.get("V2_READINGS", "").strip()
if _want:
    MODES = tuple(s.strip() for s in _want.split(",") if s.strip())
    _source = "V2_READINGS"
else:
    MODES = tuple(READINGS["in_force"])
    _source = "settings.toml, readings.in_force"
_unknown = [c for c in MODES if c not in ALL_MODES]
if _unknown:
    raise ValueError(
        f"{_source}: no such reading {_unknown}; choose from {list(ALL_MODES)}"
    )

# divisions left out of the subcortex reference of subref ('' is a structure that
# belongs to no division). The fiber tracts (fiber tracts-unassigned to cbf) and
# the ventricular system (V3 to c) are categories of the atlas, not divisions, so
# they are listed by the divisions its division term set splits them into. The two
# catch-alls go too: brain-unassigned is brain with no structure, unassigned is
# outside the brain. Left in: TH, HY, PAL, MB, P and MY
NOT_SUBCORTEX = {
    "Isocortex",
    "HPF",
    "STR",
    "OLF",
    "CTXsp",
    "CB",
    "fiber tracts-unassigned",
    "lfbs",
    "mfbs",
    "cm",
    "scwm",
    "eps",
    "cbf",
    "V3",
    "V4",
    "VL",
    "AQ",
    "c",
    "brain-unassigned",
    "unassigned",
    "",
}

# the cohorts, from per_mouse.MICE (the cohort table), each in the route's order
YOUNG_P20 = [m for m, (cohort, _, _) in MICE.items() if cohort == "young_P20"]
YOUNG_P16 = [m for m, (cohort, _, _) in MICE.items() if cohort == "young_P16"]
YOUNG_P22 = [m for m, (cohort, _, _) in MICE.items() if cohort == "young_P22"]
NAIVE = [m for m, (cohort, _, _) in MICE.items() if cohort == "naive"]
RWS = [m for m, (cohort, _, _) in MICE.items() if cohort == "rws"]
COHORTS = {
    "young": YOUNG_P20 + YOUNG_P16 + YOUNG_P22,
    "young_P20": YOUNG_P20,
    "young_P16": YOUNG_P16,
    "young_P22": YOUNG_P22,
    "naive": NAIVE,
    "rws": RWS,
    "adult": NAIVE + RWS,
}


def readings_in_force() -> str:
    """The readings in force as a run prints them, marked when V2_READINGS set them."""
    text = " ".join(MODES)
    if os.environ.get("V2_READINGS", "").strip():
        text += "  (from V2_READINGS)"
    return text


def mouse_scalars(mouse: str) -> dict[str, float]:
    """The per-brain numbers the readings divide by, cached beside the per-mouse file.

    Measured on the brain's own atlas before any warp. cortex_mean and
    subcortex_mean are plain means over their labels; z_median and z_spread
    describe the brain's distribution across structures of log2(structure mean
    / cortex mean), which is what the range match needs. Structure level, not
    voxel level, so a large structure cannot set the spread on its own; a
    structure counts from region_tables.min_vox20 voxels (settings.toml).
    """
    src = PER_MOUSE / (mouse + ".npz")
    cache = PER_MOUSE / (mouse + "_scalars.npz")

    # the cache records the date of the file it was computed from, so rerunning
    # run_per_mouse invalidates it rather than leaving a stale cortex mean
    if cache.exists():
        z = np.load(cache)
        if "src_mtime" in z.files and float(z["src_mtime"]) == src.stat().st_mtime:
            return {k: float(z[k]) for k in z.files if k != "src_mtime"}

    # structure and division of each parcellation index
    _, structure, division = structure_terms()

    # voxel count and signal sum per label, over the tissue
    ann = annotation_20(MICE[mouse][1])
    z = np.load(PER_MOUSE / (mouse + ".npz"))
    sig = z["sig"].astype(np.float32)
    tissue = z["tissue"]
    labels = ann[tissue]
    n_labels = int(ann.max()) + 1
    n = np.bincount(labels, minlength=n_labels)
    total = np.bincount(labels, weights=sig[tissue], minlength=n_labels)
    iso = [i for i in structure if division.get(i) == "Isocortex" and i < n_labels]
    sub = [
        i for i in structure if division.get(i, "") not in NOT_SUBCORTEX and i < n_labels
    ]
    cortex_mean = total[iso].sum() / n[iso].sum()
    sub_mean = total[sub].sum() / n[sub].sum()

    # labels pooled per structure name, then each structure's log2 relative level
    per_struct = {}
    for i in np.nonzero(n)[0]:
        if i == 0 or i not in structure:
            continue
        a, b = per_struct.get(structure[i], (0, 0.0))
        per_struct[structure[i]] = (a + int(n[i]), b + total[i])
    vals = np.array(
        [
            np.log2(t / c / cortex_mean)
            for c, t in per_struct.values()
            if c >= REGION_TABLES["min_vox20"] and t > 0
        ]
    )
    p10, med, p90 = np.percentile(vals, [10, 50, 90])
    out = dict(
        cortex_mean=float(cortex_mean),
        subcortex_mean=float(sub_mean),
        z_median=float(med),
        z_spread=float(max(p90 - p10, 1e-6)),
    )
    np.savez(cache, src_mtime=src.stat().st_mtime, **out)
    return out


def per_unit(num: np.ndarray, ref: np.ndarray, tissue: np.ndarray) -> np.ndarray:
    """`num` per unit of the reference channel `ref`, voxel by voxel.

    The denominator is smoothed first (readings.ref_sigma, one 20 um voxel), so a
    single dark voxel in the reference cannot blow the ratio up; the smoothing is
    normalised by `tissue`, so tissue at the edge is not divided by the black
    outside it. The ratio is clipped to +-readings.ratio_clip. Off tissue, and where
    the smoothed reference is not positive, there is no ratio: NaN, which the cohort
    volumes and the region tables (young_vs_adult.region_plot and region_groups)
    leave out of their sums and counts.
    """
    sigma = READINGS["ref_sigma"]
    ref_s = gaussian_filter(np.where(tissue, ref, 0), sigma) / np.maximum(
        gaussian_filter(tissue.astype(np.float32), sigma), 1e-3
    )
    keep = tissue & (ref_s > 0)
    r = np.where(keep, num / np.maximum(ref_s, 1e-3), np.nan)
    return np.clip(r, -READINGS["ratio_clip"], READINGS["ratio_clip"])


def finite_sums(
    labels: np.ndarray, values: np.ndarray, n_labels: int
) -> tuple[np.ndarray, np.ndarray]:
    """Sum and count of the finite `values` per label; a NaN enters neither."""
    ok = np.isfinite(values)
    sums = np.bincount(labels[ok], weights=values[ok], minlength=n_labels)
    counts = np.bincount(labels[ok], minlength=n_labels)
    return sums, counts


def mouse_modes(mouse: str) -> tuple[dict[str, np.ndarray], np.ndarray]:
    """The readings of one brain in the CCF, NaN off tissue, and its tissue mask.

    Returns ({reading: volume}, tissue) for the readings in MODES; a reading is
    computed only when it is asked for.
    """
    z = np.load(PER_MOUSE_CCF / (mouse + ".npz"))
    sig = z["sig"].astype(np.float32)
    auto = z["auto"].astype(np.float32)
    tissue = z["tissue"]

    # only sepratio needs the SEP channel, so the chain still runs on a brain
    # that run_add_sep_channel has not reached yet
    if "sepratio" in MODES and "sep" not in z.files:
        raise ValueError(
            f"{mouse}: no SEP channel in its per-mouse CCF file. Run\n"
            "  run_add_sep_channel.m for this brain, "
            "then run_per_mouse.py and run_to_ccf.py,\n"
            "  or drop the reading with V2_READINGS."
        )
    scalars = mouse_scalars(mouse)
    cref = np.where(tissue, sig / scalars["cortex_mean"], np.nan)
    subref = np.where(tissue, sig / scalars["subcortex_mean"], np.nan)
    # floor of sig / cortex mean, so log2 stays finite in the dimmest tissue
    floored = np.maximum(sig / scalars["cortex_mean"], READINGS["log2_floor"])
    zref = np.where(
        tissue,
        (np.log2(floored) - scalars["z_median"]) / scalars["z_spread"],
        np.nan,
    )

    # each reading computed only when MODES asks for it
    out = {
        "ratio": lambda: per_unit(sig, auto, tissue).astype(np.float32),
        "sepratio": lambda: per_unit(sig, z["sep"].astype(np.float32), tissue).astype(
            np.float32
        ),
        "cref": lambda: cref.astype(np.float32),
        "subref": lambda: subref.astype(np.float32),
        "zref": lambda: zref.astype(np.float32),
    }
    return ({k: out[k]() for k in MODES}, tissue)


def new_sums(volumes: dict[str, np.ndarray]) -> dict[str, list[np.ndarray]]:
    """Empty sums per reading: of the values, of their squares, and the count."""
    return {
        k: [
            np.zeros(v.shape, np.float64),
            np.zeros(v.shape, np.float64),
            np.zeros(v.shape, np.int16),
        ]
        for k, v in volumes.items()
    }


def add_mouse(acc: dict[str, list[np.ndarray]], volumes: dict[str, np.ndarray]) -> None:
    """Add one brain's readings to the sums, in place; a NaN adds nothing."""
    for k, v in volumes.items():
        w = np.nan_to_num(v)
        acc[k][0] += w
        acc[k][1] += w * w
        acc[k][2] += np.isfinite(v)


def write_stats(out: Path, acc: dict[str, list[np.ndarray]], tag: str) -> None:
    """Write <reading><tag>_{mean,sd,n}.npy, NaN where fewer than one or two mice."""
    for k, (s, ss, n) in acc.items():
        nf = np.maximum(n, 1).astype(np.float64)
        mean = np.where(n > 0, s / nf, np.nan).astype(np.float32)
        var = np.where(n > 1, (ss - s * s / nf) / np.maximum(nf - 1, 1), np.nan)
        np.save(out / f"{k}{tag}_mean.npy", mean)
        np.save(
            out / f"{k}{tag}_sd.npy",
            np.sqrt(np.maximum(var, 0)).astype(np.float32),
        )
        np.save(out / f"{k}{tag}_n.npy", n)


def main() -> None:
    """Write the mean, SD and n of every reading for every cohort, a line each.

    Twice per reading: over the whole brain, and with each brain's hemispheres
    averaged first (`_folded`, the left half).
    """
    for cohort, mice in COHORTS.items():
        out = OUT_ROOT / cohort
        out.mkdir(parents=True, exist_ok=True)

        # per reading: the sum of the values, the sum of their squares, and the count of
        # the mice with a value there (ratio and sepratio have none in tissue where the
        # smoothed reference is not positive, and a NaN must not count as a zero); the
        # same with each brain's two hemispheres averaged first
        acc = acc_folded = None
        for mouse in mice:
            modes, tissue = mouse_modes(mouse)

            # a voxel with no value on either side stays NaN, which nanmean warns about
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                folded = {k: fold(v) for k, v in modes.items()}
            if acc is None:
                acc, acc_folded = new_sums(modes), new_sums(folded)
            add_mouse(acc, modes)
            add_mouse(acc_folded, folded)

        # mean and sample SD, NaN where fewer than one or two mice have a value
        write_stats(out, acc, "")
        write_stats(out, acc_folded, "_folded")
        with open(out / "mice.txt", "w") as fh:
            fh.write("\n".join(mice) + "\n")
        n = acc["cref"][2]
        print(
            f"{cohort:10s} {len(mice):2d} mice | "
            f"voxels with n>=1 {int((n > 0).sum()):>11,d} | "
            f"with all mice {int((n == len(mice)).sum()):>11,d}",
            flush=True,
        )
