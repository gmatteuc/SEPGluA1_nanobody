"""
v2 normalisation, step 2: cohort volumes from the per-mouse files.

Two ways of expressing each mouse's background-subtracted nano, kept side by
side because they answer different questions:
  ratio   sig / auto per voxel, the autofluorescence channel as the internal
          standard (the auto denominator is smoothed by one 20 um voxel so
          single dark voxels do not blow up). Comparable across ages to the
          extent that autofluorescence is age-stable.
  cref    sig / (isocortex mean of sig), a pure scale with no intercept.
          Region ratios within a mouse are preserved exactly. Cortex-relative
          by construction, so blind to a whole-cortex change.

Per voxel the cohort gets the MEAN over the mice that have tissue there, the
SD, and N. No voxel is required to be covered by every mouse; the N map says
what each mean rests on and the downstream reader thresholds it. Both
hemispheres are kept; the hemisphere average is done at comparison time.

Output: data/comparisons_v2/<cohort>/{ratio,cref}_{mean,sd,n}.npy at 20 um
(AP, DV, ML) on the cohort's own atlas grid, plus mice.txt. Cohorts: each of
young_P20 / naive / rws, and adult = naive + rws.

  D:\\sep_histology\\code\\tools\\venv_atlas\\Scripts\\python.exe v2_cohort.py
"""

import os

import numpy as np
from scipy.ndimage import gaussian_filter

DATA = r'D:\sep_histology\data'
V2 = os.path.join(DATA, 'comparisons_v2')
PER_MOUSE = os.path.join(V2, 'per_mouse')
COHORTS = {
    'young_P20': ['MG897_SepGluA_P20', 'MG903_SepGluA_P20', 'MG913_SepGluA_P20'],
    'naive': ['CGF027_Gria1', 'CGF028_Gria1', 'CGF033_Gria1', 'CGF034_Gria1', 'CGF035_Gria1'],
    'rws': ['MG691_Gria1', 'MG692_Gria1', 'MG693_Gria1', 'MG736_Gria1', 'MG737_Gria1'],
}
COHORTS['adult'] = COHORTS['naive'] + COHORTS['rws']
RATIO_CLIP = 20.0


def mouse_modes(mouse):
    z = np.load(os.path.join(PER_MOUSE, mouse + '.npz'))
    sig = z['sig'].astype(np.float32); auto = z['auto'].astype(np.float32); tissue = z['tissue']
    auto_s = gaussian_filter(np.where(tissue, auto, 0), 1.0) / np.maximum(gaussian_filter(tissue.astype(np.float32), 1.0), 1e-3)
    ratio = np.where(tissue & (auto_s > 0), sig / np.maximum(auto_s, 1e-3), np.nan)
    ratio = np.clip(ratio, -RATIO_CLIP, RATIO_CLIP)
    cref = np.where(tissue, sig / float(z['cortex_mean']), np.nan)
    return {'ratio': ratio.astype(np.float32), 'cref': cref.astype(np.float32)}, tissue


def main():
    for cohort, mice in COHORTS.items():
        out = os.path.join(V2, cohort); os.makedirs(out, exist_ok=True)
        acc = None
        for mouse in mice:
            modes, tissue = mouse_modes(mouse)
            if acc is None:
                acc = {k: [np.zeros(v.shape, np.float64), np.zeros(v.shape, np.float64), np.zeros(v.shape, np.int16)]
                       for k, v in modes.items()}
            for k, v in modes.items():
                w = np.nan_to_num(v)
                acc[k][0] += w; acc[k][1] += w * w; acc[k][2] += tissue
            print(f'  {cohort}: {mouse} accumulated', flush=True)
        for k, (s, ss, n) in acc.items():
            nf = np.maximum(n, 1).astype(np.float64)
            mean = np.where(n > 0, s / nf, np.nan).astype(np.float32)
            var = np.where(n > 1, (ss - s * s / nf) / np.maximum(nf - 1, 1), np.nan)
            sd = np.sqrt(np.maximum(var, 0)).astype(np.float32)
            np.save(os.path.join(out, f'{k}_mean.npy'), mean)
            np.save(os.path.join(out, f'{k}_sd.npy'), sd)
            np.save(os.path.join(out, f'{k}_n.npy'), n)
        with open(os.path.join(out, 'mice.txt'), 'w') as fh:
            fh.write('\n'.join(mice) + '\n')
        print(f'{cohort}: {len(mice)} mice, voxels with n>=1: {(acc["cref"][2] > 0).sum():,}, '
              f'with all mice: {(acc["cref"][2] == len(mice)).sum():,}', flush=True)


if __name__ == '__main__':
    main()
