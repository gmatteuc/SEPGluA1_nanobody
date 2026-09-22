"""
v2, step 3: cohort volumes, built in the adult CCF from the per-mouse files.

Everything is averaged on one grid, 660 x 400 x 570 at 20 um, because v2_to_ccf
has already carried each brain there: the young through the DeMBA -> CCF
transform of its own age, the adults by placement alone. That is what lets a
pooled young group mix P20 and P16 without either age being carried by the
other's deformation field.

Two readings of the same background-subtracted signal, kept side by side:
  ratio   sig / auto per voxel (auto smoothed by one 20 um voxel so a dark
          voxel cannot blow it up): nano per unit autofluorescence, the
          internal standard. Absolute-ish across ages, to the extent that
          autofluorescence is age-stable, which is an assumption.
  cref    sig / (that mouse's isocortex mean of sig, measured natively before
          any warp): a pure scale, so region ratios within a mouse survive
          exactly. Cortex-relative by construction, hence blind to a change
          that moves the whole cortex.

Per voxel the cohort gets the mean over the mice that have tissue there, the
SD and the count. No voxel is required to be covered by every mouse; the n map
says what each mean rests on and the reader thresholds it.

Cohorts: young (P20 + P16 pooled, the group Sami wants), young_P20 (the
sensitivity check without the P16 brain), young_P16, naive, rws, and
adult = naive + rws.

Output: data/comparisons_v2/ccf/<cohort>/{ratio,cref}_{mean,sd,n}.npy + mice.txt

  D:\\sep_histology\\code\\tools\\venv_atlas\\Scripts\\python.exe v2_cohort.py
"""

import os

import numpy as np
from scipy.ndimage import gaussian_filter

from v2_per_mouse import DATA

V2 = os.path.join(DATA, 'comparisons_v2')
PER_MOUSE_CCF = os.path.join(V2, 'per_mouse_ccf')
OUT_ROOT = os.path.join(V2, 'ccf')
RATIO_CLIP = 20.0

YOUNG_P20 = ['MG897_SepGluA_P20', 'MG903_SepGluA_P20', 'MG913_SepGluA_P20',
             'MG909_SepGluA_P20', 'MG910_SepGluA_P20']
YOUNG_P16 = ['MG911_SepGluA_P16']
NAIVE = ['CGF027_Gria1', 'CGF028_Gria1', 'CGF033_Gria1', 'CGF034_Gria1', 'CGF035_Gria1']
RWS = ['MG691_Gria1', 'MG692_Gria1', 'MG693_Gria1', 'MG736_Gria1', 'MG737_Gria1']
COHORTS = {
    'young': YOUNG_P20 + YOUNG_P16,
    'young_P20': YOUNG_P20,
    'young_P16': YOUNG_P16,
    'naive': NAIVE,
    'rws': RWS,
    'adult': NAIVE + RWS,
}


def mouse_modes(mouse):
    """(dict of mode -> volume with NaN off tissue, tissue mask) for one brain in CCF."""
    z = np.load(os.path.join(PER_MOUSE_CCF, mouse + '.npz'))
    sig = z['sig'].astype(np.float32); auto = z['auto'].astype(np.float32); tissue = z['tissue']
    auto_s = gaussian_filter(np.where(tissue, auto, 0), 1.0) / np.maximum(gaussian_filter(tissue.astype(np.float32), 1.0), 1e-3)
    ratio = np.where(tissue & (auto_s > 0), sig / np.maximum(auto_s, 1e-3), np.nan)
    ratio = np.clip(ratio, -RATIO_CLIP, RATIO_CLIP)
    cref = np.where(tissue, sig / float(z['cortex_mean']), np.nan)
    return {'ratio': ratio.astype(np.float32), 'cref': cref.astype(np.float32)}, tissue


def main():
    for cohort, mice in COHORTS.items():
        out = os.path.join(OUT_ROOT, cohort); os.makedirs(out, exist_ok=True)
        acc = None
        for mouse in mice:
            modes, tissue = mouse_modes(mouse)
            if acc is None:
                acc = {k: [np.zeros(v.shape, np.float64), np.zeros(v.shape, np.float64), np.zeros(v.shape, np.int16)]
                       for k, v in modes.items()}
            for k, v in modes.items():
                w = np.nan_to_num(v)
                acc[k][0] += w; acc[k][1] += w * w; acc[k][2] += tissue
        for k, (s, ss, n) in acc.items():
            nf = np.maximum(n, 1).astype(np.float64)
            mean = np.where(n > 0, s / nf, np.nan).astype(np.float32)
            var = np.where(n > 1, (ss - s * s / nf) / np.maximum(nf - 1, 1), np.nan)
            np.save(os.path.join(out, f'{k}_mean.npy'), mean)
            np.save(os.path.join(out, f'{k}_sd.npy'), np.sqrt(np.maximum(var, 0)).astype(np.float32))
            np.save(os.path.join(out, f'{k}_n.npy'), n)
        with open(os.path.join(out, 'mice.txt'), 'w') as fh:
            fh.write('\n'.join(mice) + '\n')
        n = acc['cref'][2]
        print(f'{cohort:10s} {len(mice):2d} mice | voxels with n>=1 {int((n > 0).sum()):>11,d} | '
              f'with all mice {int((n == len(mice)).sum()):>11,d}', flush=True)


if __name__ == '__main__':
    main()
