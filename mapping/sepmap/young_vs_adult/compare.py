"""Young against adult in the adult CCF: folded maps and a per-structure table.

The cohort volumes are already on the adult CCF grid, 660 x 400 x 570 at 20 um.
This module folds the two hemispheres, compares the young cohort with the adults
reading by reading, and writes the maps, figures and tables run_compare.py lists.
An intensity reading is compared by a ratio, a signed one by a difference:

    ratio, sepratio, cref, subref:   log2(young / adult), each floored at Z_FLOOR
    zref:                            young - adult

ratio (nano over autofluorescence) is the closest to an absolute level; cref
(relative to each brain's isocortex) describes the distribution. zref is already a
position on a log scale within each brain, so its difference is the comparison.
The floor is the one volumes.cohort applies, so a map and a table agree.

The young group pools the P16, P20 and P22 brains (YOUNG); the P20 brains alone
(YOUNG_ALT) are compared too, to show what the other ages do to the answer. A
voxel is compared only where at least MIN_N_YOUNG young and MIN_N_ADULT adult
brains have tissue; the figures show the rest in grey.

Run by run_compare.py.
"""

import csv
import os
import time

import numpy as np
import nibabel as nib
from scipy.ndimage import gaussian_filter
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sepmap.config import SETTINGS
from sepmap.volumes.per_mouse import DATA, CSV_MAP
from sepmap.volumes.cohort import (
    OUT_ROOT as CCF_ROOT,
    COHORTS,
    MODES,
    SIGNED_READINGS,
    Z_FLOOR,
)

OUT = os.path.join(DATA, "comparisons_v2", "young_vs_adult")

# start of the adult registered crop along AP, in 10 um planes; added to the plane
# number in each slice title
CCF_AP0 = 180

# brains with tissue that a voxel needs to be compared, from settings.toml, the
# same keys young_vs_adult.closeup reads
MIN_N_YOUNG = SETTINGS["young_vs_adult"]["min_n_young"]
MIN_N_ADULT = SETTINGS["young_vs_adult"]["min_n_adult"]

# Gaussian sigma in 20 um voxels, applied to the log2 map only
SMOOTH = 1.0

# the pooled young cohort (P16, P20 and P22 brains), and the P20 brains alone as
# the sensitivity check
YOUNG = "young"
YOUNG_ALT = "young_P20"


def save_figure(fig, path):
    """Save `fig` as a PNG at `path` and as an EPS beside it.

    Windows refuses to overwrite a PNG that an image viewer holds open. The
    figure then goes to <name>_new.png, with a note, so a run that writes
    several figures does not lose the rest because one of them was being
    looked at.

    The EPS is what goes into a figure for a paper. PostScript has no
    transparency, so the image layers are rasterised and composited by Agg
    first; otherwise a no-data region, transparent here, would come out opaque
    black instead of showing the ground beneath it. Text, lines and axes stay
    vector, the part that has to be editable.
    """
    try:
        fig.savefig(path, dpi=105)
    except OSError:
        alt = path.replace(".png", "_new.png")
        fig.savefig(alt, dpi=105)
        print(
            f"  NOTE: {os.path.basename(path)} is open elsewhere; "
            f"wrote {os.path.basename(alt)} instead",
            flush=True,
        )

    # the EPS, with the image layers rasterised
    eps = os.path.splitext(path)[0] + ".eps"
    for ax in fig.axes:
        for im in ax.images:
            im.set_rasterized(True)
    try:
        fig.savefig(eps, dpi=105, facecolor=fig.get_facecolor(), format="eps")
    except OSError:
        print(
            f"  NOTE: {os.path.basename(eps)} is open elsewhere; "
            "the PNG was still written",
            flush=True,
        )


def fold(v):
    """Average the two hemispheres of an (AP, DV, ML) volume, ignoring NaN.

    The right half is mirrored onto the left; the result is (AP, DV, ML/2).
    """
    ml = v.shape[2]
    h = ml // 2
    left, right = v[:, :, :h], v[:, :, ml - h :][:, :, ::-1]
    return np.nanmean(np.stack([left, right]), axis=0)


def fold_n(n):
    """Fold a count map as fold does, keeping the larger count of the two sides."""
    ml = n.shape[2]
    h = ml // 2
    return np.maximum(n[:, :, :h], n[:, :, ml - h :][:, :, ::-1])


def load_cohort(cohort):
    """Folded mean maps of a cohort in the CCF, and its folded count map.

    Returns ({reading: map}, n), where n counts the brains with tissue (cref_n).
    """
    m = {k: fold(np.load(os.path.join(CCF_ROOT, cohort, f"{k}_mean.npy"))) for k in MODES}
    n = fold_n(np.load(os.path.join(CCF_ROOT, cohort, "cref_n.npy")))
    return m, n


def draw_figures(z, out):
    """Draw one slices figure per reading from the saved maps `z` into `out`.

    Also called by young_vs_adult.replot. Six coronal planes, spread over the
    planes that at least half the maximum coverage reaches; in each, the adult
    mean, the young mean and their comparison. No data (outside the atlas, or
    too few brains) is flat grey. The intensity scale stops before the white
    end of hot, so saturation reads as bright yellow and is never confused with
    empty; its range is the 99th percentile of the adult isocortex, the
    structure the question is about, so the hippocampus saturates by design.
    Signed readings use a purple-orange scale, symmetric at the 98th percentile
    of their absolute value.
    """
    from matplotlib.colors import LinearSegmentedColormap

    ann_h, both = z["annot20"], z["both"]
    inside = ann_h > 0

    # hot up to 0.82 of its range, transparent where there is no value
    hot = plt.get_cmap("hot")
    hot_cut = LinearSegmentedColormap.from_list("hot_cut", hot(np.linspace(0, 0.82, 256)))
    hot_cut.set_bad((0, 0, 0, 0))

    # zref is a position, not an intensity: its own diverging scale (purple low,
    # orange high, no green), red-blue being kept for the young-adult difference
    puor = plt.get_cmap("PuOr_r").copy()
    puor.set_bad((0, 0, 0, 0))
    grey = "#bfbfbf"

    # isocortex voxels, for the colour range
    iso_ids = set()
    with open(CSV_MAP, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            if (
                row["parcellation_term_set_name"] == "division"
                and row["parcellation_term_acronym"] == "Isocortex"
            ):
                iso_ids.add(int(row["parcellation_index"]))
    iso = both & np.isin(ann_h, list(iso_ids))

    # six planes between the first and the last that half the maximum coverage reaches
    cov = both.reshape(both.shape[0], -1).sum(1)
    ok = np.nonzero(cov > 0.5 * cov.max())[0]
    planes = [int(v) for v in np.linspace(ok[0], ok[-1], 6).round()]

    # first line of each figure's title
    what = {
        "ratio": "nano / autofluorescence, both background-subtracted",
        "sepratio": "nano / SEP, both background-subtracted  -  NOT a surface fraction, "
        "see run_sep_channel_check",
        "cref": "background-subtracted nano relative to each mouse's isocortex mean",
        "subref": "background-subtracted nano relative to the subcortex, "
        "excluding HPF and STR",
        "zref": "range-matched: position within each brain's own distribution "
        "(median 0, p90-p10 = 1), so level and dynamic range are equal across brains",
    }
    n_young = int(z["n_young_mice"])
    n_adult = int(z["n_adult_mice"])

    for m in MODES:
        adult_v, young_v, log2_v = z[f"adult_{m}"], z[f"young_{m}"], z[f"log2_{m}"]
        signed = m in SIGNED_READINGS

        # colour limits: adult isocortex for an intensity, symmetric for a signed
        # reading, and symmetric for the young-adult comparison
        vmax = (
            float(np.nanpercentile(adult_v[iso], 99))
            if not signed
            else float(np.nanpercentile(np.abs(adult_v[both]), 98))
        )
        lim2 = float(np.nanpercentile(np.abs(log2_v[both]), 98))

        # one row per plane: adult, young, comparison
        fig, axes = plt.subplots(len(planes), 3, figsize=(13.5, 3.9 * len(planes)))
        for i, zc in enumerate(planes):
            shown = both[zc]
            cmap_mean = puor if signed else hot_cut
            lim_mean = (-vmax, vmax) if signed else (0, vmax)
            diff_name = "young - adult" if signed else "log2( young / adult )"
            panels = (
                (
                    np.where(inside[zc] & np.isfinite(adult_v[zc]), adult_v[zc], np.nan),
                    f"adult (n = {n_adult})  {m}",
                    cmap_mean,
                    lim_mean,
                ),
                (
                    np.where(shown, young_v[zc], np.nan),
                    f"young (n = {n_young})  {m}",
                    cmap_mean,
                    lim_mean,
                ),
                (np.where(shown, log2_v[zc], np.nan), diff_name, "RdBu_r", (-lim2, lim2)),
            )
            for j, (im, title, cmap, lim) in enumerate(panels):
                ax = axes[i, j]

                # the atlas in flat grey under the data, so no data reads as grey
                bg = np.zeros(inside[zc].shape + (4,))
                bg[inside[zc]] = matplotlib.colors.to_rgba(grey)
                ax.imshow(bg, origin="upper", interpolation="nearest", aspect="equal")
                h = ax.imshow(
                    im,
                    cmap=cmap,
                    vmin=lim[0],
                    vmax=lim[1],
                    origin="upper",
                    interpolation="nearest",
                    aspect="equal",
                )
                ax.set_title(f"{title}   plane {zc * 2 + CCF_AP0} / 10 um", fontsize=9.5)
                ax.set_xticks([])
                ax.set_yticks([])
                plt.colorbar(
                    h, ax=ax, fraction=0.035, pad=0.01, extend="max" if j < 2 else "both"
                )

        # title: the reading, the method, and what the colour range means
        scale_note = (
            f"purple-orange scale, 0 = that brain's median structure, "
            f"+-{vmax:.2f} = its own p10-p90 spread; "
            f"red-blue on the right is the young-adult difference"
            if signed
            else f"colour range 0 to {vmax:.2f} = 99th percentile of adult isocortex; "
            "brighter is yellow, never white"
        )
        fig.suptitle(
            "\n".join(
                (
                    what[m] + ".",
                    "Hemispheres averaged; every brain carried into the adult CCF "
                    "individually.  "
                    f"Grey = no data (fewer than {MIN_N_YOUNG} young or {MIN_N_ADULT} "
                    "adults with tissue).",
                    scale_note,
                )
            ),
            fontsize=10,
        )
        fig.tight_layout(rect=(0, 0, 1, 0.975))
        save_figure(fig, os.path.join(out, f"slices_{m}.png"))
        plt.close(fig)


def main():
    """Compare young with adult and write the maps, the tables and the figures.

    The maps are compared voxel by voxel and smoothed (Gaussian, SMOOTH), with
    the voxels outside the compared mask given no weight. The tables average
    per structure, not per parcellation index: in this ontology the layers of
    an area are separate indices that share one structure name, and a table of
    areas is what anyone reads. A lookup from index to structure does that
    grouping once, and bincount then sums 18 million voxels without a Python
    loop. A group with no value at a voxel (the P20 group inside the pooled
    group's mask) is left out of its own count rather than poisoning its mean.
    Only structures with at least 100 compared voxels enter the table.
    """
    t0 = time.time()
    os.makedirs(OUT, exist_ok=True)

    # the annotation at 20 um, left half: the young volumes live on the same
    # 660-plane grid as the adults (planes 90-539 hold the adult registered crop)
    ann = np.asarray(
        nib.load(os.path.join(DATA, "atlas", "annotation_10.nii.gz")).dataobj
    )[::2, ::2, ::2]
    ann_h = ann[:, :, :285]
    inside = ann_h > 0

    # cohorts, and the voxels where both groups have enough brains
    adult, adult_n = load_cohort("adult")
    young, young_n = load_cohort(YOUNG)
    young_alt, young_alt_n = load_cohort(YOUNG_ALT)
    both = (
        inside
        & (young_n >= MIN_N_YOUNG)
        & (adult_n >= MIN_N_ADULT)
        & np.isfinite(young["cref"])
        & np.isfinite(adult["cref"])
    )
    print(
        f"voxels compared: {both.sum():,} of {inside.sum():,} inside the atlas "
        f"(young n>={MIN_N_YOUNG}: {(inside & (young_n >= MIN_N_YOUNG)).sum():,}; "
        f"adult n>={MIN_N_ADULT}: {(inside & (adult_n >= MIN_N_ADULT)).sum():,})   "
        f"{time.time() - t0:.0f} s",
        flush=True,
    )

    # young against adult per voxel, then smoothed inside the compared mask; the
    # floor is the one volumes.cohort applies, so a map and a table agree
    eps = Z_FLOOR
    log2, log2_alt = {}, {}
    for m in MODES:
        for src, dst in ((young, log2), (young_alt, log2_alt)):
            a = np.where(both, adult[m], np.nan)
            p = np.where(both, src[m], np.nan)

            # zref is already a log-scale position, so the two groups are
            # compared by difference; everything else by log2 ratio
            r = (
                (p - a)
                if m in SIGNED_READINGS
                else np.log2(np.maximum(p, eps) / np.maximum(a, eps))
            )

            # smooth r over the compared voxels only, then mask again
            w = both.astype(np.float32)
            num = gaussian_filter(np.nan_to_num(r) * w, SMOOTH)
            den = gaussian_filter(w, SMOOTH)
            dst[m] = np.where(both, num / np.maximum(den, 1e-3), np.nan).astype(
                np.float32
            )
    np.savez_compressed(
        os.path.join(OUT, "volumes_ccf20.npz"),
        annot20=ann_h,
        both=both,
        adult_n=adult_n,
        young_n=young_n,
        young_alt_n=young_alt_n,
        n_young_mice=len(COHORTS[YOUNG]),
        n_adult_mice=len(COHORTS["adult"]),
        **{f"adult_{m}": adult[m].astype(np.float32) for m in MODES},
        **{f"young_{m}": young[m].astype(np.float32) for m in MODES},
        **{f"young_alt_{m}": young_alt[m].astype(np.float32) for m in MODES},
        **{f"log2_{m}": log2[m] for m in MODES},
        **{f"log2_alt_{m}": log2_alt[m] for m in MODES},
    )

    # names, acronyms and divisions of the parcellation indices
    names, acro, divi = {}, {}, {}
    with open(CSV_MAP, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            idx = int(row["parcellation_index"])
            if row["parcellation_term_set_name"] == "structure":
                names[idx] = row["parcellation_term_name"]
                acro[idx] = row["parcellation_term_acronym"]
            elif row["parcellation_term_set_name"] == "division":
                divi[idx] = row["parcellation_term_acronym"]

    # lookup from parcellation index to structure, one entry per structure name
    struct_of = np.zeros(int(ann_h.max()) + 1, np.int64)
    struct_names = []
    seen = {}
    for idx in np.unique(ann_h):
        if idx == 0:
            continue
        nm = names.get(int(idx), f"id{idx}")
        if nm not in seen:
            seen[nm] = len(struct_names)
            struct_names.append((nm, acro.get(int(idx), ""), divi.get(int(idx), "")))
        struct_of[idx] = seen[nm]

    # mean per structure of each reading and group, over its voxels with a value
    labs = struct_of[ann_h[both].astype(np.int64)]
    n_struct = len(struct_names)
    n_vox = np.bincount(labs, minlength=n_struct)
    means = {}
    for m in MODES:
        for tag, src in (("adult", adult), ("young", young), ("young_P20", young_alt)):
            v = src[m][both]
            ok = np.isfinite(v)
            tot = np.bincount(labs[ok], weights=v[ok], minlength=n_struct)
            cnt = np.bincount(labs[ok], minlength=n_struct)
            with np.errstate(invalid="ignore", divide="ignore"):
                means[(tag, m)] = np.where(cnt > 0, tot / np.maximum(cnt, 1), np.nan)

    # one row per structure of at least 100 voxels, compared as the maps are
    rows = []
    for g in np.nonzero(n_vox >= 100)[0]:
        nm, ac, dv = struct_names[g]
        r = {"structure": nm, "acronym": ac, "division": dv, "voxels_20um": int(n_vox[g])}
        for m in MODES:
            for tag in ("adult", "young", "young_P20"):
                r[f"{tag}_{m}"] = float(means[(tag, m)][g])
            if m in SIGNED_READINGS:
                r[f"log2_{m}"] = r[f"young_{m}"] - r[f"adult_{m}"]
                r[f"log2_{m}_P20only"] = r[f"young_P20_{m}"] - r[f"adult_{m}"]
            else:
                r[f"log2_{m}"] = np.log2(
                    max(r[f"young_{m}"], eps) / max(r[f"adult_{m}"], eps)
                )
                r[f"log2_{m}_P20only"] = np.log2(
                    max(r[f"young_P20_{m}"], eps) / max(r[f"adult_{m}"], eps)
                )
        rows.append(r)
    rows.sort(key=lambda r: (r["division"], r["acronym"]))
    with open(
        os.path.join(OUT, "region_table.csv"), "w", newline="", encoding="utf-8"
    ) as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(
                {k: (f"{v:.4f}" if isinstance(v, float) else v) for k, v in r.items()}
            )

    # the cortical areas, largest cref difference first, visual and somatosensory marked
    cortex = sorted(
        [r for r in rows if r["division"] == "Isocortex"], key=lambda r: -r["log2_cref"]
    )
    lines = [f"{'area':9s} {'structure':34s} " + " ".join(f"{m:>10s}" for m in MODES)]
    for r in cortex:
        tag = (
            "  <-- visual"
            if r["acronym"].startswith("VIS")
            else ("  <-- somatosensory" if r["acronym"].startswith("SS") else "")
        )
        lines.append(
            f"{r['acronym']:9s} {r['structure'][:34]:34s} "
            + " ".join(f"{r[f'log2_{m}']:+10.2f}" for m in MODES)
            + tag
        )
    with open(os.path.join(OUT, "cortex_table.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))

    # figures, from the saved maps
    draw_figures(dict(np.load(os.path.join(OUT, "volumes_ccf20.npz"))), OUT)
    print(f"done  {time.time() - t0:.0f} s")
