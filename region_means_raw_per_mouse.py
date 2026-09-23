"""
SUPERSEDED by v2_per_mouse.py + v2_region_plot.py (Sep 2026).
Same idea -- per-mouse region means from the raw stacks, no P6bis -- but it
takes tissue from P6bis background masks, which call a partial section all
background, and it knows only the three P20 brains of the first round.
v2 takes tissue from the autofluorescence channel instead and handles any
age. Kept as the record of how the raw-stack route was validated.

Per-mouse, per-region mean of the RAW nano channel, for the young and adult
cohorts, on each cohort's own atlas.

Why raw. P6bis maps every mouse onto its cohort's cortex with a slope AND an
intercept. Ratios between regions do not survive an affine map, so any
"region A relative to reference B" read off P6bis output is distorted, and
distorted differently for each mouse. The P2bis-equalised volumes that P5
stacks (nano_4d<tag>.mat) have only had per-slice scaling within each mouse,
which preserves region ratios. So reduce those to a table -- one row per
mouse and region -- and let the choice of reference be a division applied
afterwards, with per-mouse spread for a proper test.

Why no transform. A region-wise comparison only needs the label of every
voxel, and the DeMBA labels were remapped to Allen parcellation_index when
that atlas was built, so P20 (DeMBA grid) and adults (CCF grid) share the
ontology. Each cohort is measured on its own atlas; nothing is warped.

Tissue: P6bis's background mask AND raw > 0 (the registered volume is zero
where a mouse's sections have run out, and the background detector calls an
empty plane all-tissue). Regions with fewer than MIN_VOX tissue voxels in a
mouse are left blank for that mouse rather than estimated from a sliver.

Output: data/comparisons/young_P20_vs_adult_nano/region_means_raw_per_mouse.csv
  cohort, mouse, parcellation_index, structure, acronym, division, n_vox, mean_raw

  D:\\sep_histology\\code\\tools\\venv_atlas\\Scripts\\python.exe region_means_raw_per_mouse.py
"""

import csv
import os
import time

import h5py
import nibabel as nib
import numpy as np

DATA = r'D:\sep_histology\data'
OUT = os.path.join(DATA, 'comparisons', 'young_P20_vs_adult_nano')
CSV_MAP = os.path.join(DATA, 'atlas', 'parcellation_to_parcellation_term_membership.csv')
MIN_VOX = 2000     # 10 um voxels: 2 nl

CHANNEL = os.environ.get('SEP_CHANNEL', 'nano')
assert CHANNEL in ('nano', 'auto'), CHANNEL
YOUNG = ['MG897_SepGluA_P20', 'MG903_SepGluA_P20', 'MG913_SepGluA_P20']
COHORTS = [
    # cohort, raw 4D file (None -> per-mouse registered tiff), variable, bkg-mask file, mice, atlas
    ('young_P20', os.path.join(DATA, 'young', 'nano_4d_P20.mat') if CHANNEL == 'nano' else None, 'nano_4d',
     os.path.join(DATA, 'young', 'nano_4d_normalized_bkgmask_P20.mat'), YOUNG, 'demba_p20'),
    ('naive', os.path.join(DATA, 'naive', CHANNEL + '_4d.mat'), CHANNEL + '_4d',
     os.path.join(DATA, 'naive', 'nano_4d_normalized_bkgmask.mat'),
     ['CGF027_Gria1', 'CGF028_Gria1', 'CGF033_Gria1', 'CGF034_Gria1', 'CGF035_Gria1'], 'ccf'),
    ('rws', os.path.join(DATA, 'rws', CHANNEL + '_4d.mat'), CHANNEL + '_4d',
     os.path.join(DATA, 'rws', 'nano_4d_normalized_bkgmask.mat'),
     ['MG691_Gria1', 'MG692_Gria1', 'MG693_Gria1', 'MG736_Gria1', 'MG737_Gria1'], 'ccf'),
]
TIFF = {'nano': 'chan02_NANO.tiff', 'auto': 'chan03_AUTO.tiff'}


class TiffPages:
    """Per-mouse registered tiff as an (AP, DV)-page reader: page i is ML plane i."""
    def __init__(self, mouse):
        import tifffile
        self.t = tifffile.TiffFile(os.path.join(DATA, 'young', mouse, 'lightsuite', 'volume_registered', TIFF[CHANNEL]))
        self.n = len(self.t.pages)

    def page(self, i):
        return np.asarray(self.t.pages[i].asarray(), dtype=np.float32)      # (AP, DV) already


def labels_on_registered_grid(atlas_key):
    """Annotation on the (AP, DV, ML) 10 um-equivalent grid of the registered volumes."""
    if atlas_key == 'ccf':
        ann = np.asarray(nib.load(os.path.join(DATA, 'atlas', 'annotation_10.nii.gz')).dataobj)
        return ann[179:1079]                                             # [180 1079], 900 planes
    ann = np.asarray(nib.load(os.path.join(DATA, 'atlas_demba_p20', 'annotation_10.nii.gz')).dataobj)
    ann = ann[62:559]                                                    # [63 559], 497 planes at 20 um
    return np.repeat(np.repeat(np.repeat(ann, 2, 0), 2, 1), 2, 2)        # x2 nearest -> 994 x 800 x 1140


def main():
    t0 = time.time()
    names, acro, divi = {}, {}, {}
    with open(CSV_MAP, newline='', encoding='utf-8') as fh:
        for row in csv.DictReader(fh):
            idx = int(row['parcellation_index'])
            if row['parcellation_term_set_name'] == 'structure':
                names[idx] = row['parcellation_term_name']; acro[idx] = row['parcellation_term_acronym']
            elif row['parcellation_term_set_name'] == 'division':
                divi[idx] = row['parcellation_term_acronym']

    rows = []
    label_cache = {}
    for cohort, raw_path, var, bk_path, mice, atlas_key in COHORTS:
        if atlas_key not in label_cache:
            label_cache[atlas_key] = labels_on_registered_grid(atlas_key)
        lab = label_cache[atlas_key]                                     # (AP, DV, ML)
        nlab = int(lab.max()) + 1
        fb = h5py.File(bk_path, 'r')
        B = fb['recomputed_bkg_mask_4d']                                   # (mouse, ML, DV, AP)
        if raw_path is None:
            fr, R = None, None
        else:
            fr = h5py.File(raw_path, 'r'); R = fr[var]                        # (mouse, ML, DV, AP)
            assert R.shape[1:] == lab.shape[::-1], (R.shape, lab.shape)
        for m, mouse in enumerate(mice):
            src = TiffPages(mouse) if R is None else None
            sums = np.zeros(nlab, np.float64); cnts = np.zeros(nlab, np.int64)
            bg_samples = []
            for i in range(lab.shape[2]):                                  # one ML page at a time
                v = src.page(i) if src else np.asarray(R[m, i], dtype=np.float32).T   # (AP, DV)
                bgm = np.asarray(B[m, i]).T
                tissue = (bgm == 0) & (v > 0)
                l = lab[:, :, i][tissue]
                sums += np.bincount(l, weights=v[tissue], minlength=nlab)
                cnts += np.bincount(l, minlength=nlab)
                if i % 4 == 0:
                    bg_samples.append(v[(bgm == 1) & (v > 0)][::4])
            bg = np.concatenate(bg_samples)
            bg_median = float(np.median(bg)) if bg.size else float('nan')
            rows.append((cohort, mouse, -1, 'BACKGROUND (off-tissue median)', 'bg', 'bg', int(bg.size * 16), bg_median))
            for idx in np.nonzero(cnts >= MIN_VOX)[0]:
                if idx == 0:
                    continue
                rows.append((cohort, mouse, int(idx), names.get(int(idx), f'id{idx}'), acro.get(int(idx), ''),
                             divi.get(int(idx), ''), int(cnts[idx]), sums[idx] / cnts[idx]))
            print(f'{cohort:10s} {mouse:20s} {int((cnts >= MIN_VOX).sum()):4d} regions   bg median {bg_median:.0f}   {time.time() - t0:.0f} s', flush=True)
        fb.close()
        if fr is not None:
            fr.close()

    out = os.path.join(OUT, 'region_means_raw_per_mouse' + ('' if CHANNEL == 'nano' else '_' + CHANNEL) + '.csv')
    with open(out, 'w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(['cohort', 'mouse', 'parcellation_index', 'structure', 'acronym', 'division', 'n_vox', 'mean_raw'])
        w.writerows(rows)
    print(f'wrote {out}: {len(rows)} rows, {time.time() - t0:.0f} s')


if __name__ == '__main__':
    main()
