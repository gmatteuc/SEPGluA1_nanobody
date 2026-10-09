"""The coronal plane the figures of the ISH line draw, and what goes on it.

The guided figures show one CCF plane, ish_figures.plane in 10 um planes (700:
cortex, hippocampus and thalamus in one section). Three grids meet there:

    labels       the adult CCF annotation on the 20 um grid of the adult crop
                 (volumes.per_mouse.annotation_20), which starts at 10 um plane
                 CROP_START, so 10 um plane k is crop plane (k - CROP_START) / 2
    nano map     the cohort mean of cref (each adult over its own isocortex mean)
                 on the full 20 um CCF grid (volumes.cohort), plane k / 2, shown
                 where at least videos.min_n.adult adults have tissue
    ISH section  an Allen 200 um grid, section k / 20, each 20 um pixel taking its
                 nearest 200 um voxel, so nothing is interpolated; sections set
                 missing by the section QC stay missing
    one adult    a raw channel of one adult's per-mouse file (counts above the
                 channel's off-tissue background), on the adult crop like the
                 labels, NaN outside that brain's tissue

Nothing here computes a result: these are the images a figure puts beside the
structure values it shows (comparison_rows, the genes of figure 05), and the atlas
of the QC sheets (sheet_atlas).

Used by the run scripts of the ISH line that draw coronal panels.
"""

import nibabel as nib
import numpy as np
import pandas as pd

from sepmap import structures
from sepmap.config import DATA, SETTINGS
from sepmap.ish import section_qc
from sepmap.ish.regions import read_energy
from sepmap.volumes.cohort import OUT_ROOT
from sepmap.volumes.per_mouse import OUT as PER_MOUSE
from sepmap.volumes.per_mouse import annotation_20, atlas_grid, structure_terms

# the adults with tissue a voxel of the nano map needs to be shown, as in the adult
# cohort's videos
VIDEOS = SETTINGS["videos"]

# the first 10 um CCF plane of the adult crop
CROP_START = atlas_grid("ccf")[1][0]

ADULT_MAP = OUT_ROOT / "adult"


def crop_index(plane: int) -> int:
    """The plane of the adult crop's 20 um grid that holds a 10 um CCF plane."""
    return (plane - CROP_START) // 2


def label_plane(plane: int) -> np.ndarray:
    """CCF parcellation indices of one 10 um plane, on the 20 um grid, (DV, ML)."""
    return np.asarray(annotation_20("ccf")[crop_index(plane)])


def channel_planes(
    mouse: str, plane: int, channels: tuple[str, ...] = ("sig", "sep", "auto")
) -> dict[str, np.ndarray]:
    """One adult's raw channels on one 10 um plane, {channel: (DV, ML)}.

    Counts above each channel's background, NaN outside the brain's tissue mask.
    """
    z = np.load(PER_MOUSE / f"{mouse}.npz")
    k = crop_index(plane)
    tissue = np.asarray(z["tissue"][k])
    return {
        c: np.where(tissue, np.asarray(z[c][k], dtype=np.float32), np.nan)
        for c in channels
    }


def nano_plane(plane: int) -> np.ndarray:
    """The adult cohort's cref on one plane, (DV, ML); NaN where few adults reach."""
    mean = np.load(ADULT_MAP / "cref_mean.npy", mmap_mode="r")[plane // 2]
    n = np.load(ADULT_MAP / "cref_n.npy", mmap_mode="r")[plane // 2]
    shown = n >= VIDEOS["min_n"]["adult"]
    return np.where(shown, np.asarray(mean, dtype=np.float32), np.nan)


def ish_plane(
    experiment_id: str,
    plane_of_section: str,
    flagged: list[int],
    plane: int,
    shape: tuple[int, int],
) -> np.ndarray:
    """One Allen experiment's energy at a 10 um plane, on the 20 um grid, (DV, ML).

    `plane_of_section` (coronal or sagittal) says along which axis the sections in
    `flagged`, those the QC set missing, lie; `shape` is the label plane's. The 200
    um grid samples the atlas every 20th 10 um plane, so each 20 um pixel takes
    its nearest 200 um voxel.
    """
    vol = read_energy(experiment_id)
    axis = section_qc.SECTION_AXIS[plane_of_section]
    vol = section_qc.apply_flags(vol, flagged, axis)
    section = vol[plane // 20]
    rows = np.clip(np.rint(np.arange(shape[0]) / 10).astype(int), 0, section.shape[0] - 1)
    cols = np.clip(np.rint(np.arange(shape[1]) / 10).astype(int), 0, section.shape[1] - 1)
    return section[np.ix_(rows, cols)]


def sheet_atlas() -> tuple[np.ndarray, np.ndarray]:
    """The template and the structures of the whole CCF on the 20 um grid, (AP, DV, ML).

    Ten times finer than the ISH grid, so a QC sheet draws the borders of structures
    over a section; the structures as codes (structures.name_volume), the layers of
    an area one label.
    """
    template = nib.load(DATA / "atlas" / "average_template_10.nii.gz")
    template = np.asarray(template.dataobj)[::2, ::2, ::2]
    ann = np.asarray(nib.load(DATA / "atlas" / "annotation_10.nii.gz").dataobj)
    names, _, _ = structure_terms()
    labels, _ = structures.name_volume(ann[::2, ::2, ::2], names)
    return template, labels


def comparison_rows(
    nano: pd.DataFrame,
    merged: dict[str, dict[str, float]],
    table: pd.DataFrame,
    lab: np.ndarray,
    plane: int,
    genes: tuple[str, ...],
) -> list[dict]:
    """The genes of figure 05: each one's ISH section on the plane, beside its rows.

    `nano` is the nano map's rows of gene_ranking.csv by gene, `merged` the merged
    profiles and `table` the gene table. The section shown is P9's own coronal
    experiment when it is used, else the gene's first coronal one; its caption says
    which, and how many the gene has.
    """
    rows = []
    for gene in genes:
        used = table[(table["symbol"] == gene) & ~table["excluded"]]
        coronal = used[used["plane"] == "coronal"]
        own = coronal[coronal["p9_experiment"]]
        shown = (own if len(own) else coronal).iloc[0]
        flagged = [int(k) for k in str(shown["flagged_sections"]).split()]
        image = ish_plane(shown["experiment_id"], "coronal", flagged, plane, lab.shape)
        caption = (
            f"Allen experiment {shown['experiment_id']}, one of the gene's "
            f"{len(used)} (200 um)"
        )
        if len(used) == 1:
            caption = f"Allen experiment {shown['experiment_id']}, its only one (200 um)"
        rows.append(
            dict(
                symbol=gene,
                image=image,
                caption=caption,
                profile=merged[gene],
                ranking=nano.loc[gene],
                n_all=len(nano),
            )
        )
    return rows
