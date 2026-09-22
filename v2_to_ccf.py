"""
v2, step 2: put every brain on the adult CCF grid, one brain at a time.

Why per mouse rather than per cohort mean. The young group now spans two ages,
P20 and P16, whose atlases have different AP extents (497 and 478 planes of the
705-plane DeMBA canvas), so their volumes cannot be averaged before they are
brought to a common space. Warping each brain individually also means the
cohort mean, its spread and its n map are all computed on the same grid, and a
pooled group can mix ages without any of them being carried by another age's
deformation field.

  young   block-averaged registered volume -> placed back into its age's full
          DeMBA canvas -> brainglobe_ccf_translator from age_PND to P56 in
          allen_mouse -> 660 x 400 x 570 at 20 um. About 2.5 min per volume,
          three volumes per mouse (sig, auto, tissue).
  adults  already registered to the CCF crop [180 1079] at 10 um, i.e. exactly
          planes 90..539 of the same 20 um CCF grid, so they are only placed,
          never warped. Warping them would blur them for nothing.

The tissue mask travels as a mask (nearest neighbour) and is re-thresholded
after the transform, so a warped voxel is tissue only if it came from tissue.

Output: data/comparisons_v2/per_mouse_ccf/<mouse>.npz with sig, auto (float16)
and tissue (bool) on the 660 x 400 x 570 CCF grid at 20 um, plus the scalars
the per-mouse file carried (backgrounds, cortex mean, cohort, age).

  D:\\sep_histology\\code\\tools\\venv_atlas\\Scripts\\python.exe v2_to_ccf.py [mouse ...]
"""

import os
import re
import sys
import time

import numpy as np

from v2_per_mouse import DATA, MICE, OUT as PER_MOUSE, atlas_grid

OUT = os.path.join(DATA, 'comparisons_v2', 'per_mouse_ccf')
DEMBA_SHAPE = (705, 400, 570)        # every DeMBA age shares this canvas at 20 um
CCF_SHAPE = (660, 400, 570)          # allen_mouse at 20 um
CCF_CROP = (90, 540)                 # the adult registered crop, 10 um planes 180..1079


def to_ccf(vol_demba_full, age, is_mask=False):
    """(AP, DV, ML) on the full DeMBA canvas at 20 um -> (AP, DV, ML) in the adult CCF.

    demba_dev_mouse wants (ML, DV, AP); transform() works in place and hands
    back (AP, DV, ML). This is the call validated on the annotation.
    """
    from brainglobe_ccf_translator import Volume
    v = Volume(values=np.ascontiguousarray(np.transpose(vol_demba_full, (2, 1, 0))),
               space='demba_dev_mouse', voxel_size_micron=20, age_PND=age,
               segmentation_file=is_mask)
    v.transform(target_age=56, target_space='allen_mouse')
    return np.asarray(v.values, dtype=np.float32)


def main(mice):
    os.makedirs(OUT, exist_ok=True)
    for mouse in mice:
        t0 = time.time()
        cohort, atlas_key = MICE[mouse][:2]
        z = np.load(os.path.join(PER_MOUSE, mouse + '.npz'))
        sig = z['sig'].astype(np.float32); auto = z['auto'].astype(np.float32); tissue = z['tissue']
        lo, hi = atlas_grid(atlas_key)[1]

        if atlas_key == 'ccf':
            # no warp: drop the adult crop into the full CCF grid at 20 um
            out = {}
            for name, arr, fill in (('sig', sig, 0.0), ('auto', auto, 0.0), ('tissue', tissue, False)):
                full = np.full(CCF_SHAPE, fill, arr.dtype)
                full[CCF_CROP[0]:CCF_CROP[1]] = arr
                out[name] = full
            age = 56
        else:
            age = int(re.search(r'p(\d+)$', atlas_key).group(1))
            out = {}
            for name, arr, is_mask in (('sig', sig, False), ('auto', auto, False), ('tissue', tissue, True)):
                full = np.zeros(DEMBA_SHAPE, np.float32)
                full[lo - 1:hi] = np.where(tissue, arr, 0) if name != 'tissue' else tissue.astype(np.float32)
                out[name] = to_ccf(full, age, is_mask=is_mask)
                print(f'  {mouse} {name} warped P{age} -> CCF   {time.time() - t0:.0f} s', flush=True)
            out['tissue'] = out['tissue'] > 0.5
            # a warped value only counts where the warped mask says tissue
            for name in ('sig', 'auto'):
                out[name] = np.where(out['tissue'], out[name], 0.0)

        np.savez_compressed(os.path.join(OUT, mouse + '.npz'),
                            sig=out['sig'].astype(np.float16), auto=out['auto'].astype(np.float16),
                            tissue=out['tissue'], cohort=cohort, atlas=atlas_key, age=age,
                            bg_nano=z['bg_nano'], bg_auto=z['bg_auto'], cortex_mean=z['cortex_mean'])
        print(f'{mouse:20s} {cohort:10s} P{age:<3d} tissue {int(out["tissue"].sum()):>10,d} voxels in CCF   '
              f'{time.time() - t0:.0f} s', flush=True)


if __name__ == '__main__':
    main(sys.argv[1:] or list(MICE))
