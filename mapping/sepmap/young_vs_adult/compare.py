"""Young against adult in the adult CCF: folded maps and a per-structure table.

The cohort volumes are already on the adult CCF grid, 660 x 400 x 570 at 20 um.
This module folds the two hemispheres, compares the young cohort with the adults
reading by reading, and writes the maps, figures and tables run_compare.py lists.
An intensity reading is compared by a ratio, a signed one by a difference:

    ratio, sepratio, cref, subref:   log2(young / adult), each floored at
                                     readings.log2_floor (settings.toml)
    zref:                            young - adult

ratio (nano over autofluorescence) is the closest to an absolute level; cref
(relative to each brain's isocortex) describes the distribution. zref is already a
position on a log scale within each brain, so its difference is the comparison.
The floor is the one volumes.cohort applies, so a map and a table agree.

The young group pools the P16, P20 and P22 brains (YOUNG); the P20 brains alone
(YOUNG_ALT) are compared too, to show what the other ages do to the answer. A
voxel is compared only where at least young_vs_adult.min_n_young young and
min_n_adult adult brains have tissue; the figures show the rest in grey.

Run by run_compare.py.
"""

import csv
import time
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
from matplotlib.colors import Colormap
from scipy.ndimage import gaussian_filter

from sepmap.config import SETTINGS
from sepmap.plotting import NO_DATA_GREY, hot_cut, save_figure, transparent_bad
from sepmap.volumes.cohort import COHORTS, MODES, SIGNED_READINGS
from sepmap.volumes.cohort import OUT_ROOT as CCF_ROOT
from sepmap.volumes.per_mouse import DATA, isocortex_ids, structure_terms
from sepmap.young_vs_adult.hemispheres import fold, fold_count

READINGS = SETTINGS["readings"]
YOUNG_VS_ADULT = SETTINGS["young_vs_adult"]

OUT = DATA / "comparisons_v2" / "young_vs_adult"

# the pooled young cohort (P16, P20 and P22 brains), and the P20 brains alone as
# the sensitivity check
YOUNG = "young"
YOUNG_ALT = "young_P20"

# first line of each slices figure's title
READING_TITLES = {
    "ratio": "nano / autofluorescence, both background-subtracted",
    "sepratio": "nano / SEP, both background-subtracted  -  NOT a surface fraction, "
    "see run_sep_channel_check",
    "cref": "background-subtracted nano relative to each mouse's isocortex mean",
    "subref": "background-subtracted nano relative to the subcortex, "
    "excluding HPF and STR",
    "zref": "range-matched: position within each brain's own distribution "
    "(median 0, p90-p10 = 1), so level and dynamic range are equal across brains",
}


def load_cohort(cohort: str) -> tuple[dict[str, np.ndarray], np.ndarray]:
    """Folded mean maps of a cohort in the CCF, and its folded count map.

    Returns ({reading: map}, n), where n counts the brains with tissue (cref_n).
    """
    maps = {
        reading: fold(np.load(CCF_ROOT / cohort / f"{reading}_mean.npy"))
        for reading in MODES
    }
    n = fold_count(np.load(CCF_ROOT / cohort / "cref_n.npy"))
    return maps, n


def slice_planes(compared: np.ndarray) -> list[int]:
    """Six planes from the first to the last that half the maximum coverage reaches."""
    coverage = compared.reshape(compared.shape[0], -1).sum(1)
    covered = np.nonzero(coverage > 0.5 * coverage.max())[0]
    planes = [int(v) for v in np.linspace(covered[0], covered[-1], 6).round()]
    return planes


def draw_plane(
    axes: np.ndarray, i: int, zc: int, panels: tuple, inside: np.ndarray
) -> None:
    """Draw row `i` of a slices figure, plane `zc`: (image, title, colormap, limits)."""
    for j, (im, title, cmap, lim) in enumerate(panels):
        ax = axes[i, j]

        # the atlas in flat grey under the data, so no data reads as grey
        bg = np.zeros(inside[zc].shape + (4,))
        bg[inside[zc]] = matplotlib.colors.to_rgba(NO_DATA_GREY)
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
        ax.set_title(f"{title}   CCF plane {zc * 2} / 10 um", fontsize=9.5)
        ax.set_xticks([])
        ax.set_yticks([])
        plt.colorbar(
            h, ax=ax, fraction=0.035, pad=0.01, extend="max" if j < 2 else "both"
        )


def slices_title(reading: str, signed: bool, vmax: float) -> str:
    """The title of a slices figure: the reading, the method, what the colours mean."""
    min_n_young = YOUNG_VS_ADULT["min_n_young"]
    min_n_adult = YOUNG_VS_ADULT["min_n_adult"]
    if signed:
        scale_note = (
            "purple-orange scale, 0 = that brain's median structure, "
            f"+-{vmax:.2f} = its own p10-p90 spread; "
            "red-blue on the right is the young-adult difference"
        )
    else:
        scale_note = (
            f"colour range 0 to {vmax:.2f} = 99th percentile of adult isocortex; "
            "brighter is yellow, never white"
        )
    return "\n".join(
        (
            READING_TITLES[reading] + ".",
            "Hemispheres averaged; every brain carried into the adult CCF "
            "individually.  "
            f"Grey = no data (fewer than {min_n_young} young or {min_n_adult} "
            "adults with tissue).",
            scale_note,
        )
    )


def draw_reading(
    reading: str,
    maps: dict[str, np.ndarray],
    planes: list[int],
    inside: np.ndarray,
    compared: np.ndarray,
    isocortex: np.ndarray,
    hot: Colormap,
    puor: Colormap,
    n_young: int,
    n_adult: int,
    out: Path,
) -> None:
    """Draw and save slices_<reading>.png: a row per plane, adult, young, comparison."""
    adult_v = maps[f"adult_{reading}"]
    young_v = maps[f"young_{reading}"]
    log2_v = maps[f"log2_{reading}"]
    signed = reading in SIGNED_READINGS

    # colour limits: adult isocortex for an intensity, symmetric for a signed
    # reading, and symmetric for the young-adult comparison
    if not signed:
        vmax = float(np.nanpercentile(adult_v[isocortex], 99))
    else:
        vmax = float(np.nanpercentile(np.abs(adult_v[compared]), 98))
    diff_lim = float(np.nanpercentile(np.abs(log2_v[compared]), 98))

    # one row per plane: adult, young, comparison
    fig, axes = plt.subplots(len(planes), 3, figsize=(13.5, 3.9 * len(planes)))
    for i, zc in enumerate(planes):
        shown = compared[zc]
        cmap_mean = puor if signed else hot
        lim_mean = (-vmax, vmax) if signed else (0, vmax)
        diff_name = "young - adult" if signed else "log2( young / adult )"
        panels = (
            (
                np.where(inside[zc] & np.isfinite(adult_v[zc]), adult_v[zc], np.nan),
                f"adult (n = {n_adult})  {reading}",
                cmap_mean,
                lim_mean,
            ),
            (
                np.where(shown, young_v[zc], np.nan),
                f"young (n = {n_young})  {reading}",
                cmap_mean,
                lim_mean,
            ),
            (
                np.where(shown, log2_v[zc], np.nan),
                diff_name,
                "RdBu_r",
                (-diff_lim, diff_lim),
            ),
        )
        draw_plane(axes, i, zc, panels, inside)

    # title: the reading, the method, and what the colour range means
    fig.suptitle(slices_title(reading, signed, vmax), fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.975))
    save_figure(fig, out / f"slices_{reading}.png", dpi=105)
    plt.close(fig)


def draw_figures(maps: dict[str, np.ndarray], out: Path) -> None:
    """Draw one slices figure per reading from the saved `maps` into `out`.

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
    annotation_left, compared = maps["annot20"], maps["both"]
    inside = annotation_left > 0

    # hot up to 0.82 of its range, transparent where there is no value
    hot = hot_cut()

    # zref is a position, not an intensity: its own diverging scale (purple low,
    # orange high, no green), red-blue being kept for the young-adult difference
    puor = transparent_bad("PuOr_r")

    # isocortex voxels, for the colour range
    iso_ids = isocortex_ids()
    isocortex = compared & np.isin(annotation_left, list(iso_ids))

    # the planes drawn, and the brains of each group
    planes = slice_planes(compared)
    n_young = int(maps["n_young_mice"])
    n_adult = int(maps["n_adult_mice"])

    for reading in MODES:
        draw_reading(
            reading,
            maps,
            planes,
            inside,
            compared,
            isocortex,
            hot,
            puor,
            n_young,
            n_adult,
            out,
        )


def compared_voxels(
    inside: np.ndarray,
    young: dict[str, np.ndarray],
    young_n: np.ndarray,
    adult: dict[str, np.ndarray],
    adult_n: np.ndarray,
    t0: float,
) -> np.ndarray:
    """The voxels where both groups have enough brains with a value; prints the counts."""
    min_n_young = YOUNG_VS_ADULT["min_n_young"]
    min_n_adult = YOUNG_VS_ADULT["min_n_adult"]
    compared = (
        inside
        & (young_n >= min_n_young)
        & (adult_n >= min_n_adult)
        & np.isfinite(young["cref"])
        & np.isfinite(adult["cref"])
    )
    print(
        f"voxels compared: {compared.sum():,} of {inside.sum():,} inside the atlas "
        f"(young n>={min_n_young}: {(inside & (young_n >= min_n_young)).sum():,}; "
        f"adult n>={min_n_adult}: {(inside & (adult_n >= min_n_adult)).sum():,})   "
        f"{time.time() - t0:.0f} s",
        flush=True,
    )
    return compared


def contrast_maps(
    compared: np.ndarray,
    adult: dict[str, np.ndarray],
    young: dict[str, np.ndarray],
    young_alt: dict[str, np.ndarray],
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    """Young against adult per voxel, smoothed inside the compared mask, per reading.

    Returns the maps of the pooled young group and of the P20 brains alone. The
    floor is the one volumes.cohort applies, so a map and a table agree.
    """
    floor = READINGS["log2_floor"]
    log2, log2_alt = {}, {}
    for reading in MODES:
        for src, dst in ((young, log2), (young_alt, log2_alt)):
            adult_v = np.where(compared, adult[reading], np.nan)
            young_v = np.where(compared, src[reading], np.nan)

            # zref is already a log-scale position, so the two groups are
            # compared by difference; everything else by log2 ratio
            if reading in SIGNED_READINGS:
                contrast = young_v - adult_v
            else:
                contrast = np.log2(
                    np.maximum(young_v, floor) / np.maximum(adult_v, floor)
                )

            # smooth the contrast over the compared voxels only, then mask again
            weights = compared.astype(np.float32)
            num = gaussian_filter(
                np.nan_to_num(contrast) * weights, YOUNG_VS_ADULT["smooth"]
            )
            den = gaussian_filter(weights, YOUNG_VS_ADULT["smooth"])
            dst[reading] = np.where(compared, num / np.maximum(den, 1e-3), np.nan).astype(
                np.float32
            )
    return log2, log2_alt


def write_maps(
    annotation_left: np.ndarray,
    compared: np.ndarray,
    adult: dict[str, np.ndarray],
    adult_n: np.ndarray,
    young: dict[str, np.ndarray],
    young_n: np.ndarray,
    young_alt: dict[str, np.ndarray],
    young_alt_n: np.ndarray,
    log2: dict[str, np.ndarray],
    log2_alt: dict[str, np.ndarray],
) -> None:
    """Write volumes_ccf20.npz: the folded maps, the comparisons and the counts."""
    np.savez_compressed(
        OUT / "volumes_ccf20.npz",
        annot20=annotation_left,
        both=compared,
        adult_n=adult_n,
        young_n=young_n,
        young_alt_n=young_alt_n,
        n_young_mice=len(COHORTS[YOUNG]),
        n_adult_mice=len(COHORTS["adult"]),
        **{f"adult_{reading}": adult[reading].astype(np.float32) for reading in MODES},
        **{f"young_{reading}": young[reading].astype(np.float32) for reading in MODES},
        **{
            f"young_alt_{reading}": young_alt[reading].astype(np.float32)
            for reading in MODES
        },
        **{f"log2_{reading}": log2[reading] for reading in MODES},
        **{f"log2_alt_{reading}": log2_alt[reading] for reading in MODES},
    )


def structure_means(
    annotation_left: np.ndarray,
    compared: np.ndarray,
    adult: dict[str, np.ndarray],
    young: dict[str, np.ndarray],
    young_alt: dict[str, np.ndarray],
) -> tuple[list[tuple[str, str, str]], np.ndarray, dict[tuple[str, str], np.ndarray]]:
    """Each structure's compared voxels and mean of each reading and group.

    Returns the structures as (name, acronym, division), their voxel counts, and
    {(group, reading): mean per structure}, NaN where a group has no value.
    """
    # names, acronyms and divisions of the parcellation indices
    names, acro, divi = structure_terms()

    # lookup from parcellation index to structure, one entry per structure name
    struct_of = np.zeros(int(annotation_left.max()) + 1, np.int64)
    struct_names = []
    seen = {}
    for idx in np.unique(annotation_left):
        if idx == 0:
            continue
        name = names.get(int(idx), f"id{idx}")
        if name not in seen:
            seen[name] = len(struct_names)
            struct_names.append((name, acro.get(int(idx), ""), divi.get(int(idx), "")))
        struct_of[idx] = seen[name]

    # mean per structure of each reading and group, over its voxels with a value
    structure_of_voxel = struct_of[annotation_left[compared].astype(np.int64)]
    n_struct = len(struct_names)
    n_vox = np.bincount(structure_of_voxel, minlength=n_struct)
    means = {}
    for reading in MODES:
        for tag, src in (("adult", adult), ("young", young), ("young_P20", young_alt)):
            v = src[reading][compared]
            ok = np.isfinite(v)
            total = np.bincount(structure_of_voxel[ok], weights=v[ok], minlength=n_struct)
            count = np.bincount(structure_of_voxel[ok], minlength=n_struct)
            with np.errstate(invalid="ignore", divide="ignore"):
                means[(tag, reading)] = np.where(
                    count > 0, total / np.maximum(count, 1), np.nan
                )
    return struct_names, n_vox, means


def region_rows(
    struct_names: list[tuple[str, str, str]],
    n_vox: np.ndarray,
    means: dict[tuple[str, str], np.ndarray],
) -> list[dict]:
    """One row per structure with enough voxels, compared as the maps are, sorted."""
    floor = READINGS["log2_floor"]
    rows = []
    for structure in np.nonzero(n_vox >= YOUNG_VS_ADULT["min_table_vox20"])[0]:
        name, acronym, division = struct_names[structure]
        row = {
            "structure": name,
            "acronym": acronym,
            "division": division,
            "voxels_20um": int(n_vox[structure]),
        }
        for reading in MODES:
            for tag in ("adult", "young", "young_P20"):
                row[f"{tag}_{reading}"] = float(means[(tag, reading)][structure])
            if reading in SIGNED_READINGS:
                row[f"log2_{reading}"] = row[f"young_{reading}"] - row[f"adult_{reading}"]
                row[f"log2_{reading}_P20only"] = (
                    row[f"young_P20_{reading}"] - row[f"adult_{reading}"]
                )
            else:
                row[f"log2_{reading}"] = np.log2(
                    max(row[f"young_{reading}"], floor)
                    / max(row[f"adult_{reading}"], floor)
                )
                row[f"log2_{reading}_P20only"] = np.log2(
                    max(row[f"young_P20_{reading}"], floor)
                    / max(row[f"adult_{reading}"], floor)
                )
        rows.append(row)
    rows.sort(key=lambda row: (row["division"], row["acronym"]))
    return rows


def write_region_table(rows: list[dict]) -> None:
    """Write region_table.csv, floats to four decimals."""
    with open(OUT / "region_table.csv", "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {k: (f"{v:.4f}" if isinstance(v, float) else v) for k, v in row.items()}
            )


def cortex_table(rows: list[dict]) -> None:
    """Write and print the cortical areas, largest cref difference first.

    Visual and somatosensory areas are marked; written to cortex_table.txt.
    """
    cortex = sorted(
        [row for row in rows if row["division"] == "Isocortex"],
        key=lambda row: -row["log2_cref"],
    )
    lines = [
        f"{'area':9s} {'structure':34s} "
        + " ".join(f"{reading:>10s}" for reading in MODES)
    ]
    for row in cortex:
        if row["acronym"].startswith("VIS"):
            tag = "  <-- visual"
        elif row["acronym"].startswith("SS"):
            tag = "  <-- somatosensory"
        else:
            tag = ""
        lines.append(
            f"{row['acronym']:9s} {row['structure'][:34]:34s} "
            + " ".join(f"{row[f'log2_{reading}']:+10.2f}" for reading in MODES)
            + tag
        )
    with open(OUT / "cortex_table.txt", "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))


def main() -> None:
    """Compare young with adult and write the maps, the tables and the figures.

    The maps are compared voxel by voxel and smoothed (Gaussian,
    young_vs_adult.smooth), with the voxels outside the compared mask given no
    weight. The tables average per structure, not per parcellation index: in this
    ontology the layers of an area are separate indices that share one structure
    name, and a table of areas is what anyone reads. A lookup from index to
    structure does that grouping once, and bincount then sums 18 million voxels
    without a Python loop. A group with no value at a voxel (the P20 group inside
    the pooled group's mask) is left out of its own count rather than poisoning its
    mean. Only structures with at least young_vs_adult.min_table_vox20 compared
    voxels enter the table.
    """
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)

    # the annotation at 20 um, left half: the young volumes live on the same
    # 660-plane grid as the adults (planes 90-539 hold the adult registered crop)
    annotation = np.asarray(nib.load(DATA / "atlas" / "annotation_10.nii.gz").dataobj)[
        ::2, ::2, ::2
    ]
    annotation_left = annotation[:, :, :285]
    inside = annotation_left > 0

    # cohorts, and the voxels where both groups have enough brains
    adult, adult_n = load_cohort("adult")
    young, young_n = load_cohort(YOUNG)
    young_alt, young_alt_n = load_cohort(YOUNG_ALT)
    compared = compared_voxels(inside, young, young_n, adult, adult_n, t0)

    # young against adult per voxel, smoothed inside the compared mask, and the maps
    log2, log2_alt = contrast_maps(compared, adult, young, young_alt)
    write_maps(
        annotation_left,
        compared,
        adult,
        adult_n,
        young,
        young_n,
        young_alt,
        young_alt_n,
        log2,
        log2_alt,
    )

    # the table per structure, and the cortical areas printed
    struct_names, n_vox, means = structure_means(
        annotation_left, compared, adult, young, young_alt
    )
    rows = region_rows(struct_names, n_vox, means)
    write_region_table(rows)
    cortex_table(rows)

    # figures, from the saved maps
    draw_figures(dict(np.load(OUT / "volumes_ccf20.npz")), OUT)
    print(f"done  {time.time() - t0:.0f} s")
