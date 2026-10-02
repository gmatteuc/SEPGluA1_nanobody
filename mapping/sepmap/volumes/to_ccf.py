"""Every brain onto the adult CCF grid at 20 um, one brain at a time.

Each brain is carried on its own rather than as part of a cohort mean: the
young brains come at several ages whose atlases have different AP extents (497
planes of the 705-plane DeMBA canvas at P20, 478 at P16), so their volumes
cannot be averaged before they are brought to a common space. Carrying each
brain also means the cohort mean, its spread and its n map are all computed on
the same grid, and a pooled group can mix ages without any of them being
carried by another age's deformation field.

  young   block-averaged registered volume -> placed back into its age's full
          DeMBA canvas -> brainglobe_ccf_translator from age_PND to P56 in
          allen_mouse -> 660 x 400 x 570 at 20 um. About 2.5 min per volume,
          four volumes per mouse where SEP is there (sig, auto, sep, tissue).
  adults  already registered to the CCF crop [180 1079] at 10 um, which is
          exactly planes 90..539 of the same 20 um CCF grid, so they are only
          placed, never warped. Warping them would blur them for nothing.

The tissue mask travels as a mask (nearest neighbour) and is re-thresholded
after the transform, so a warped voxel is tissue only if it came from tissue.

Output: data/comparisons_v2/per_mouse_ccf/<mouse>.npz with sig, auto and sep
(float16) and tissue (bool) on the 660 x 400 x 570 CCF grid at 20 um, plus the
scalars the per-mouse file carried (backgrounds, cortex mean, cohort, age).
SEP rides exactly the channels it will be divided into, through the same
transform in the same call, so nothing can drift between them.

Run by run_to_ccf.py.
"""

import os
import re
import time

import numpy as np

from sepmap.volumes.per_mouse import DATA, MICE, atlas_grid
from sepmap.volumes.per_mouse import OUT as PER_MOUSE

OUT = os.path.join(DATA, "comparisons_v2", "per_mouse_ccf")

# every DeMBA age shares this canvas at 20 um
DEMBA_SHAPE = (705, 400, 570)

# allen_mouse at 20 um
CCF_SHAPE = (660, 400, 570)

# the adult registered crop, 10 um planes 180..1079
CCF_CROP = (90, 540)


def to_ccf(vol_demba_full: np.ndarray, age: int, is_mask: bool = False) -> np.ndarray:
    """Carry an (AP, DV, ML) volume on the full DeMBA canvas into the adult CCF.

    `vol_demba_full` is at 20 um on the canvas of age `age` (postnatal days);
    the result is (AP, DV, ML) in allen_mouse at P56, float32. A mask
    (`is_mask`) is carried by nearest neighbour. demba_dev_mouse wants (ML, DV,
    AP); transform() works in place and hands back (AP, DV, ML). This is the
    call validated on the annotation.
    """
    # imported here: only the young brains need it, and venv_flat does not have it
    from brainglobe_ccf_translator import Volume

    v = Volume(
        values=np.ascontiguousarray(np.transpose(vol_demba_full, (2, 1, 0))),
        space="demba_dev_mouse",
        voxel_size_micron=20,
        age_PND=age,
        segmentation_file=is_mask,
    )
    v.transform(target_age=56, target_space="allen_mouse")
    return np.asarray(v.values, dtype=np.float32)


def main(mice: list[str]) -> None:
    """Write the CCF file of each of `mice`, with one printed line each."""
    os.makedirs(OUT, exist_ok=True)
    for mouse in mice:
        t0 = time.time()
        cohort, atlas_key = MICE[mouse][:2]
        z = np.load(os.path.join(PER_MOUSE, mouse + ".npz"))
        sig = z["sig"].astype(np.float32)
        auto = z["auto"].astype(np.float32)
        tissue = z["tissue"]
        chans = [("sig", sig), ("auto", auto)]
        if "sep" in z.files:
            chans.append(("sep", z["sep"].astype(np.float32)))
        lo, hi = atlas_grid(atlas_key)[1]

        if atlas_key == "ccf":
            # no warp: drop the adult crop into the full CCF grid at 20 um
            out = {}
            volumes = [(n, a, 0.0) for n, a in chans] + [("tissue", tissue, False)]
            for name, arr, fill in volumes:
                full = np.full(CCF_SHAPE, fill, arr.dtype)
                full[CCF_CROP[0] : CCF_CROP[1]] = arr
                out[name] = full
            age = 56
        else:
            # the young: back into the full DeMBA canvas of their age, then warped
            age = int(re.search(r"p(\d+)$", atlas_key).group(1))
            out = {}
            volumes = [(n, a, False) for n, a in chans] + [("tissue", tissue, True)]
            for name, arr, is_mask in volumes:
                full = np.zeros(DEMBA_SHAPE, np.float32)
                if name != "tissue":
                    full[lo - 1 : hi] = np.where(tissue, arr, 0)
                else:
                    full[lo - 1 : hi] = tissue.astype(np.float32)
                out[name] = to_ccf(full, age, is_mask=is_mask)
                print(
                    f"  {mouse} {name} warped P{age} -> CCF   {time.time() - t0:.0f} s",
                    flush=True,
                )
            out["tissue"] = out["tissue"] > 0.5

            # a warped value only counts where the warped mask says tissue
            for name, _ in chans:
                out[name] = np.where(out["tissue"], out[name], 0.0)

        np.savez_compressed(
            os.path.join(OUT, mouse + ".npz"),
            tissue=out["tissue"],
            cohort=cohort,
            atlas=atlas_key,
            age=age,
            bg_nano=z["bg_nano"],
            bg_auto=z["bg_auto"],
            cortex_mean=z["cortex_mean"],
            **{n: out[n].astype(np.float16) for n, _ in chans},
        )
        print(
            f"{mouse:20s} {cohort:10s} P{age:<3d} "
            f"tissue {int(out['tissue'].sum()):>10,d} voxels in CCF   "
            f"{time.time() - t0:.0f} s",
            flush=True,
        )
