"""
v2 normalisation, step 1: per mouse, from the raw registered stacks.

What was wrong with the P6bis -> P8 route for a cross-age question, in order
of damage:
  - the tissue mask came from the nano intensity (select_background_pixels),
    which flags a whole plane as background when a section covers only part
    of it, and as tissue when the plane is empty;
  - the per-mouse map was affine (slope + intercept), which does not preserve
    ratios between regions;
  - abs() of the normalised values, and cohort means over voxels where some
    mice had no tissue at all.

What this does instead, per mouse, at 20 um (2x2x2 block mean of the
registered 10 um-equivalent grid):
  tissue   auto channel above its off-tissue level (median + 4 MAD), AND the
           registered nano is non-zero, AND inside the atlas brain. The nano
           intensity never enters the mask.
  sig      nano minus the mouse's scalar off-tissue background (raw counts)
  auto     auto minus its own off-tissue background
  plus two scalars: the isocortex mean of sig (the pure-scale cortex
  reference for the 'cref' mode downstream) and the two backgrounds.

Off-tissue samples are taken from planes where the atlas says brain but the
registered nano is zero on less than half of it (i.e. planes a section
actually reached), outside the atlas brain mask, so a degenerate bg mask
cannot pollute them.

Output: data/comparisons_v2/per_mouse/<mouse>.npz  (sig, auto: float16;
tissue: bool; scalars). Nothing under data/comparisons is touched.

  D:\\sep_histology\\code\\tools\\venv_atlas\\Scripts\\python.exe v2_per_mouse.py [mouse ...]
"""

import csv
import os
import sys
import time

import h5py
import nibabel as nib
import numpy as np

DATA = r'D:\sep_histology\data'
OUT = os.path.join(DATA, 'comparisons_v2', 'per_mouse')
CSV_MAP = os.path.join(DATA, 'atlas', 'parcellation_to_parcellation_term_membership.csv')
MAD_K = 4.0

MICE = {  # mouse -> (cohort, atlas, nano source, auto source, index in 4D stack)
    **{m: ('young_P20', 'demba_p20', 'tiff', 'tiff', None)
       for m in ('MG897_SepGluA_P20', 'MG903_SepGluA_P20', 'MG913_SepGluA_P20')},
    **{m: ('naive', 'ccf', os.path.join(DATA, 'naive', 'nano_4d.mat'), os.path.join(DATA, 'naive', 'auto_4d.mat'), i)
       for i, m in enumerate(('CGF027_Gria1', 'CGF028_Gria1', 'CGF033_Gria1', 'CGF034_Gria1', 'CGF035_Gria1'))},
    **{m: ('rws', 'ccf', os.path.join(DATA, 'rws', 'nano_4d.mat'), os.path.join(DATA, 'rws', 'auto_4d.mat'), i)
       for i, m in enumerate(('MG691_Gria1', 'MG692_Gria1', 'MG693_Gria1', 'MG736_Gria1', 'MG737_Gria1'))},
}


def annotation_20(atlas_key):
    """Labels on the 20 um (AP, DV, ML) grid that the 2x2x2 block mean of the registered stack lands on."""
    if atlas_key == 'ccf':
        ann = np.asarray(nib.load(os.path.join(DATA, 'atlas', 'annotation_10.nii.gz')).dataobj)[179:1079]
        return ann[::2, ::2, ::2]                                        # 450 x 400 x 570
    ann = np.asarray(nib.load(os.path.join(DATA, 'atlas_demba_p20', 'annotation_10.nii.gz')).dataobj)
    return ann[62:559]                                                   # 497 x 400 x 570


class Source:
    """Yields (AP, DV) pages for ML index i, from a 4D .mat stack or a registered tiff."""
    def __init__(self, path, idx, mouse, chan):
        if path == 'tiff':
            import tifffile
            name = {'nano': 'chan02_NANO.tiff', 'auto': 'chan03_AUTO.tiff'}[chan]
            self.t = tifffile.TiffFile(os.path.join(DATA, 'young', mouse, 'lightsuite', 'volume_registered', name))
            self.n = len(self.t.pages); self.h = None
        else:
            self.f = h5py.File(path, 'r'); self.h = self.f[os.path.basename(path)[:-4]]; self.idx = idx
            self.n = self.h.shape[1]; self.t = None

    def page(self, i):
        if self.t is not None:
            return np.asarray(self.t.pages[i].asarray(), dtype=np.float32)
        return np.asarray(self.h[self.idx, i], dtype=np.float32).T


def block2(a):
    """2x2 block mean of a 2D page with even dims."""
    return a.reshape(a.shape[0] // 2, 2, a.shape[1] // 2, 2).mean(axis=(1, 3))


def main(mice):
    iso = set()
    with open(CSV_MAP, newline='', encoding='utf-8') as fh:
        for row in csv.DictReader(fh):
            if row['parcellation_term_set_name'] == 'division' and row['parcellation_term_acronym'] == 'Isocortex':
                iso.add(int(row['parcellation_index']))
    os.makedirs(OUT, exist_ok=True)
    anns = {}
    for mouse in mice:
        t0 = time.time()
        cohort, atlas_key, pn, pa, idx = MICE[mouse]
        if atlas_key not in anns:
            anns[atlas_key] = annotation_20(atlas_key)
        ann = anns[atlas_key]; brain = ann > 0                            # (AP, DV, ML) at 20 um
        nano, auto = Source(pn, idx, mouse, 'nano'), Source(pa, idx, mouse, 'auto')
        n_ap, n_dv, n_ml = ann.shape
        sig = np.zeros((n_ap, n_dv, n_ml), np.float32)
        aut = np.zeros((n_ap, n_dv, n_ml), np.float32)
        nz = np.zeros((n_ap, n_dv, n_ml), np.float32)                     # fraction of the 8 fine voxels with nano != 0
        for j in range(n_ml):
            s = a = z = None
            for i in (2 * j, 2 * j + 1):
                v = nano.page(i); u = auto.page(i)
                s = block2(v) if s is None else s + block2(v)
                a = block2(u) if a is None else a + block2(u)
                z = block2((v != 0).astype(np.float32)) if z is None else z + block2((v != 0).astype(np.float32))
            sig[:, :, j] = s / 2; aut[:, :, j] = a / 2; nz[:, :, j] = z / 2
        # planes a section reached: nano non-zero on at least half the atlas brain
        reached = np.array([(nz[k][brain[k]] > 0).mean() > 0.5 if brain[k].any() else False for k in range(n_ap)])
        off = (~brain) & (nz > 0) & reached[:, None, None]                # off-tissue but imaged
        bg_n = float(np.median(sig[off])); bg_a = float(np.median(aut[off]))
        mad_a = float(np.median(np.abs(aut[off] - bg_a))) * 1.4826
        tissue = brain & (nz >= 0.5) & (aut > bg_a + MAD_K * mad_a) & reached[:, None, None]
        sig -= bg_n; aut -= bg_a
        cortex = tissue & np.isin(ann, list(iso))
        cortex_mean = float(sig[cortex].mean())
        np.savez_compressed(os.path.join(OUT, mouse + '.npz'),
                            sig=sig.astype(np.float16), auto=aut.astype(np.float16), tissue=tissue,
                            reached=reached, bg_nano=bg_n, bg_auto=bg_a, mad_auto=mad_a,
                            cortex_mean=cortex_mean, cohort=cohort, atlas=atlas_key)
        print(f'{mouse:20s} {cohort:10s} bg nano {bg_n:6.0f}  bg auto {bg_a:5.0f} (mad {mad_a:4.0f})  '
              f'tissue {100 * tissue.sum() / brain.sum():5.1f}% of atlas brain, planes reached {reached.sum()}/{n_ap}  '
              f'cortex mean {cortex_mean:6.0f}   {time.time() - t0:.0f} s', flush=True)


if __name__ == '__main__':
    main(sys.argv[1:] or list(MICE))
