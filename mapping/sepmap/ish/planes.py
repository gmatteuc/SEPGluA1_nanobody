"""The coronal plane the figures of the ISH line draw, and what goes on it.

The guided figures show one CCF plane, ish_figures.plane in 10 um planes (700:
cortex, hippocampus and thalamus in one section). Three grids meet there:

    labels       the adult CCF annotation on the 20 um grid of the adult crop
                 (volumes.per_mouse.annotation_20), which starts at 10 um plane
                 CROP_START, so 10 um plane k is crop plane (k - CROP_START) / 2
    nano map     the cohort mean of cref (each adult over its own isocortex mean)
                 on the full 20 um CCF grid (volumes.cohort), plane k / 2, shown
                 where at least MIN_BRAINS adults have tissue
    ISH section  an Allen 200 um grid, section k / 20, each 20 um pixel taking its
                 nearest 200 um voxel, so nothing is interpolated; sections set
                 missing by the section QC stay missing

Nothing here computes a result: these are the images a figure puts beside the
structure values it shows.

Used by the run scripts of the ISH line that draw coronal panels.
"""

import numpy as np

from sepmap.config import DATA
from sepmap.ish import section_qc
from sepmap.ish.regions import read_energy
from sepmap.volumes.per_mouse import annotation_20

# the first 10 um CCF plane of the adult crop (volumes.per_mouse.atlas_grid)
CROP_START = 180

# adults with tissue a voxel of the nano map needs to be shown: half the ten
MIN_BRAINS = 5

ADULT_MAP = DATA / "comparisons_v2" / "ccf" / "adult"


def label_plane(plane: int) -> np.ndarray:
    """CCF parcellation indices of one 10 um plane, on the 20 um grid, (DV, ML)."""
    return np.asarray(annotation_20("ccf")[(plane - CROP_START) // 2])


def nano_plane(plane: int) -> np.ndarray:
    """The adult cohort's cref on one plane, (DV, ML); NaN where few adults reach."""
    mean = np.load(ADULT_MAP / "cref_mean.npy", mmap_mode="r")[plane // 2]
    n = np.load(ADULT_MAP / "cref_n.npy", mmap_mode="r")[plane // 2]
    return np.where(n >= MIN_BRAINS, np.asarray(mean, dtype=np.float32), np.nan)


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
