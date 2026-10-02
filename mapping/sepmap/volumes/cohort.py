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
            expressed. The green channel in this fixed, cleared tissue is
            mostly autofluorescence (rho 0.79 +- 0.04 against the auto channel
            in all ten adults; dynamic range 0.95 log2 against auto's 1.07 and
            nano's 1.93; adult.sep_channel_check), so this reading tracks ratio
            at rho 0.89 to 0.97 within every mouse. It is not a surface fraction.
  cref      sig / that mouse's isocortex mean of sig, measured before any warp:
            a pure scale, so region ratios within a mouse survive exactly; blind
            to a change that moves the whole cortex.
  subref    sig / that mouse's subcortex mean, without the divisions in
            NOT_SUBCORTEX (cortex, HPF, STR and others): a reference that is
            neither the cortex nor the two structures that dominate the scale.
  zref      log2(sig / cortex mean), minus that brain's median over structures,
            divided by its own p90-p10 spread, so level and dynamic range are
            the same in every brain; the pup brain is flatter (spread 0.89 +-
            0.25 log2 against 1.79 +- 0.24). The only reading not linear in the
            signal: its maps show a position within a range, and the compression
            it removes may itself be the finding, so it is read beside cref.

Per voxel a cohort gets the mean over the mice with tissue there, the SD and the
count. No voxel needs every mouse; the n map says what each mean rests on, and
the reader thresholds it. Cohorts: young (P16, P20 and P22 pooled, the main
comparison), young_P20 (the P20 brains alone, the sensitivity check), young_P16,
young_P22, naive, rws, and adult (naive and rws).

Output: data/comparisons_v2/ccf/<cohort>/<reading>_{mean,sd,n}.npy and mice.txt.

Run by run_cohort.py; its cohorts and readings are imported across the package.
"""

import csv
import os

import numpy as np
from scipy.ndimage import gaussian_filter

from sepmap.config import SETTINGS
from sepmap.volumes.per_mouse import CSV_MAP, DATA, MICE, annotation_20
from sepmap.volumes.per_mouse import OUT as PER_MOUSE

V2 = os.path.join(DATA, "comparisons_v2")
PER_MOUSE_CCF = os.path.join(V2, "per_mouse_ccf")
OUT_ROOT = os.path.join(V2, "ccf")

# ratio and sepratio clipped to +-20, so a near-empty reference voxel cannot
# dominate a mean
RATIO_CLIP = 20.0

# floor of sig / cortex mean, so log2 stays finite in the dimmest tissue
Z_FLOOR = SETTINGS["readings"]["log2_floor"]
MODES = ("ratio", "sepratio", "cref", "subref", "zref")

# the readings that are a signed position rather than an intensity: young and adult
# are compared by difference, not log2 ratio, on a diverging scale centred on zero.
# Defined once, in settings.toml, because one script getting it wrong would look like
# a result, not a bug; the chain imports it from here, and young_vs_adult.closeup
# reads the same key.
SIGNED_READINGS = tuple(SETTINGS["readings"]["signed"])

# V2_READINGS, a comma-separated subset of MODES ($env:V2_READINGS = 'ratio,cref' in
# PowerShell), drops the other readings from every figure, table and video without
# touching the data; young_vs_adult.region_plot and region_groups follow it too.
_want = os.environ.get("V2_READINGS", "").strip()
if _want:
    chosen = tuple(s.strip() for s in _want.split(",") if s.strip())
    unknown = [c for c in chosen if c not in MODES]
    if unknown:
        raise SystemExit(
            f"V2_READINGS: no such reading {unknown}; choose from {list(MODES)}"
        )
    MODES = chosen

# divisions left out of the subcortex reference of subref ('' is a structure that
# belongs to no division)
NOT_SUBCORTEX = {
    "Isocortex",
    "HPF",
    "STR",
    "OLF",
    "CTXsp",
    "fiber tracts",
    "VS",
    "CB",
    "",
}

# the cohorts, from the same mice as per_mouse.MICE
YOUNG_P20 = [
    "MG897_SepGluA_P20",
    "MG903_SepGluA_P20",
    "MG913_SepGluA_P20",
    "MG909_SepGluA_P20",
    "MG910_SepGluA_P20",
]
YOUNG_P16 = ["MG911_SepGluA_P16"]
YOUNG_P22 = ["MG904_SepGluA_P22"]
NAIVE = ["CGF027_Gria1", "CGF028_Gria1", "CGF033_Gria1", "CGF034_Gria1", "CGF035_Gria1"]
RWS = ["MG691_Gria1", "MG692_Gria1", "MG693_Gria1", "MG736_Gria1", "MG737_Gria1"]
COHORTS = {
    "young": YOUNG_P20 + YOUNG_P16 + YOUNG_P22,
    "young_P20": YOUNG_P20,
    "young_P16": YOUNG_P16,
    "young_P22": YOUNG_P22,
    "naive": NAIVE,
    "rws": RWS,
    "adult": NAIVE + RWS,
}


def mouse_scalars(mouse: str) -> dict[str, float]:
    """The per-brain numbers the readings divide by, cached beside the per-mouse file.

    Measured on the brain's own atlas before any warp. cortex_mean and
    subcortex_mean are plain means over their labels; z_median and z_spread
    describe the brain's distribution across structures of log2(structure mean
    / cortex mean), which is what the range match needs. Structure level, not
    voxel level, so a large structure cannot set the spread on its own; a
    structure counts from 250 voxels.
    """
    src = os.path.join(PER_MOUSE, mouse + ".npz")
    cache = os.path.join(PER_MOUSE, mouse + "_scalars.npz")

    # the cache records the date of the file it was computed from, so rerunning
    # run_per_mouse invalidates it rather than leaving a stale cortex mean
    if os.path.exists(cache):
        z = np.load(cache)
        if "src_mtime" in z.files and float(z["src_mtime"]) == os.path.getmtime(src):
            return {k: float(z[k]) for k in z.files if k != "src_mtime"}

    # structure and division of each parcellation index
    stru, divi = {}, {}
    with open(CSV_MAP, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            i = int(r["parcellation_index"])
            if r["parcellation_term_set_name"] == "structure":
                stru[i] = r["parcellation_term_acronym"]
            elif r["parcellation_term_set_name"] == "division":
                divi[i] = r["parcellation_term_acronym"]

    # voxel count and signal sum per label, over the tissue
    ann = annotation_20(MICE[mouse][1])
    z = np.load(os.path.join(PER_MOUSE, mouse + ".npz"))
    sig = z["sig"].astype(np.float32)
    tissue = z["tissue"]
    lab = ann[tissue]
    nlab = int(ann.max()) + 1
    n = np.bincount(lab, minlength=nlab)
    tot = np.bincount(lab, weights=sig[tissue], minlength=nlab)
    iso = [i for i in stru if divi.get(i) == "Isocortex" and i < nlab]
    sub = [i for i in stru if divi.get(i, "") not in NOT_SUBCORTEX and i < nlab]
    cortex_mean = tot[iso].sum() / n[iso].sum()
    sub_mean = tot[sub].sum() / n[sub].sum()

    # labels pooled per structure name, then each structure's log2 relative level
    per_struct = {}
    for i in np.nonzero(n)[0]:
        if i == 0 or i not in stru:
            continue
        a, b = per_struct.get(stru[i], (0, 0.0))
        per_struct[stru[i]] = (a + int(n[i]), b + tot[i])
    vals = np.array(
        [
            np.log2(t / c / cortex_mean)
            for c, t in per_struct.values()
            if c >= 250 and t > 0
        ]
    )
    p10, med, p90 = np.percentile(vals, [10, 50, 90])
    out = dict(
        cortex_mean=float(cortex_mean),
        subcortex_mean=float(sub_mean),
        z_median=float(med),
        z_spread=float(max(p90 - p10, 1e-6)),
    )
    np.savez(cache, src_mtime=os.path.getmtime(src), **out)
    return out


def mouse_modes(mouse: str) -> tuple[dict[str, np.ndarray], np.ndarray]:
    """The readings of one brain in the CCF, NaN off tissue, and its tissue mask.

    Returns ({reading: volume}, tissue) for the readings in MODES; a reading is
    computed only when it is asked for.
    """
    z = np.load(os.path.join(PER_MOUSE_CCF, mouse + ".npz"))
    sig = z["sig"].astype(np.float32)
    auto = z["auto"].astype(np.float32)
    tissue = z["tissue"]

    # only sepratio needs the SEP channel, so the chain still runs on a brain
    # that run_add_sep_channel has not reached yet
    if "sepratio" in MODES and "sep" not in z.files:
        raise SystemExit(
            f"{mouse}: no SEP channel in its per-mouse CCF file. Run\n"
            f"  run_add_sep_channel.m for this brain, "
            f"then run_per_mouse.py and run_to_ccf.py,\n"
            f"  or drop the reading with V2_READINGS."
        )
    sc = mouse_scalars(mouse)

    def per_unit(ref):
        """Divide sig by a reference channel, voxel by voxel.

        The denominator is smoothed by one 20 um voxel first, so a single dark
        voxel in the reference cannot blow the ratio up; the smoothing is
        normalised by the mask, so tissue at the edge is not divided by the
        black outside it.
        """
        ref_s = gaussian_filter(np.where(tissue, ref, 0), 1.0) / np.maximum(
            gaussian_filter(tissue.astype(np.float32), 1.0), 1e-3
        )
        r = np.where(tissue & (ref_s > 0), sig / np.maximum(ref_s, 1e-3), np.nan)
        return np.clip(r, -RATIO_CLIP, RATIO_CLIP).astype(np.float32)

    cref = np.where(tissue, sig / sc["cortex_mean"], np.nan)
    subref = np.where(tissue, sig / sc["subcortex_mean"], np.nan)
    floored = np.maximum(sig / sc["cortex_mean"], Z_FLOOR)
    zref = np.where(tissue, (np.log2(floored) - sc["z_median"]) / sc["z_spread"], np.nan)

    # each reading computed only when MODES asks for it
    out = {
        "ratio": lambda: per_unit(auto),
        "sepratio": lambda: per_unit(z["sep"].astype(np.float32)),
        "cref": lambda: cref.astype(np.float32),
        "subref": lambda: subref.astype(np.float32),
        "zref": lambda: zref.astype(np.float32),
    }
    return ({k: out[k]() for k in MODES}, tissue)


def main() -> None:
    """Write the mean, SD and n of every reading for every cohort, a line each."""
    for cohort, mice in COHORTS.items():
        out = os.path.join(OUT_ROOT, cohort)
        os.makedirs(out, exist_ok=True)

        # per reading: the sum of the values, the sum of their squares, the count
        acc = None
        for mouse in mice:
            modes, tissue = mouse_modes(mouse)
            if acc is None:
                acc = {
                    k: [
                        np.zeros(v.shape, np.float64),
                        np.zeros(v.shape, np.float64),
                        np.zeros(v.shape, np.int16),
                    ]
                    for k, v in modes.items()
                }
            for k, v in modes.items():
                w = np.nan_to_num(v)
                acc[k][0] += w
                acc[k][1] += w * w
                acc[k][2] += tissue

        # mean and sample SD, NaN where fewer than one or two mice have tissue
        for k, (s, ss, n) in acc.items():
            nf = np.maximum(n, 1).astype(np.float64)
            mean = np.where(n > 0, s / nf, np.nan).astype(np.float32)
            var = np.where(n > 1, (ss - s * s / nf) / np.maximum(nf - 1, 1), np.nan)
            np.save(os.path.join(out, f"{k}_mean.npy"), mean)
            np.save(
                os.path.join(out, f"{k}_sd.npy"),
                np.sqrt(np.maximum(var, 0)).astype(np.float32),
            )
            np.save(os.path.join(out, f"{k}_n.npy"), n)
        with open(os.path.join(out, "mice.txt"), "w") as fh:
            fh.write("\n".join(mice) + "\n")
        n = acc["cref"][2]
        print(
            f"{cohort:10s} {len(mice):2d} mice | "
            f"voxels with n>=1 {int((n > 0).sum()):>11,d} | "
            f"with all mice {int((n == len(mice)).sum()):>11,d}",
            flush=True,
        )
