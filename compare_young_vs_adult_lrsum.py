"""
SUPERSEDED by v2_to_ccf.py + v2_cohort.py + v2_compare.py (Sep 2026).
It reads the P8 mean_lr_sum volumes, which carry the P6bis affine map and
the abs() of it, so its region ratios are distorted; and it warps a cohort
MEAN into the CCF, which cannot work once the young group spans two ages.
Kept because the coverage masks and the DeMBA -> CCF call were first
worked out here, and the numbers in the early notes came from it.

Rapid P20 vs adult comparison of the hemisphere-folded nano LR-sum.

Takes the cohort mean LR-sum volumes P8 writes (mean_lr_sum_*.mat), carries the
P20 one from DeMBA space into the adult CCF with the CCF Translator
deformation fields (the transform validated on the annotation: 0.945 voxel
agreement), and compares the two on the same 20 um grid.

The two cohorts were each normalised within themselves (P6bis, then P8's
common factor), so their absolute units are not comparable. Each volume is
therefore divided by its own brain-wide mean inside the analysis range before
comparing, and the comparison is a LOG2 RATIO of those relative levels: 0 =
same share of the brain's signal, +1 = twice the share in P20, -1 = half.
That is a statement about the spatial DISTRIBUTION of the label, not about
absolute expression.

Outputs, in data/comparisons/young_P20_vs_adult_nano/:
  region_log2ratio.csv       per structure: adult rel, P20 rel, log2 ratio, voxels
  slices_adult_vs_P20.png    three coronal planes: adult | P20 in CCF | log2 ratio
  volumes_ccf20.npz          both relative volumes + masks on the CCF 20 um half grid
  P20_lrsum_in_ccf20.nii.gz  the transformed P20 mean LR-sum, full CCF 20 um grid

Run with the atlas environment:
  D:\\sep_histology\\code\\tools\\venv_atlas\\Scripts\\python.exe compare_young_vs_adult_lrsum.py
"""

import csv
import os
import time
from collections import defaultdict

import h5py
import nibabel as nib
import numpy as np
from scipy.ndimage import gaussian_filter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from brainglobe_ccf_translator import Volume

DATA = r'D:\sep_histology\data'
P20_MAT = os.path.join(DATA, 'comparisons', 'merged_young_P20_nano', 'mean_lr_sum_merged_young_P20_nano_nosmooth.mat')
ADU_MAT = os.path.join(DATA, 'comparisons', 'merged_naive_rws_nano', 'mean_lr_sum_merged_naive_rws_nano_nosmooth.mat')
CCF_ANNOT = os.path.join(DATA, 'atlas', 'annotation_10.nii.gz')
CSV_MAP = os.path.join(DATA, 'atlas', 'parcellation_to_parcellation_term_membership.csv')
OUT = os.path.join(DATA, 'comparisons', 'young_P20_vs_adult_nano')
os.makedirs(OUT, exist_ok=True)

DEMBA_APLIMS = (63, 559)      # 1-based, 20 um, as in get_atlas('demba_p20')
CCF_APLIMS = (180, 1079)      # 1-based, 10 um, as in get_atlas('ccf')
ANALYSIS_10UM = (100, 700)    # adult planes used for region means in P8
SMOOTH_SIGMA_20UM = 1.5       # light smoothing before the ratio, in 20 um voxels


def load_mat(path):
    """mean_lr_sum and masks as (AP, DV, ML_half) float32 -- MATLAB order."""
    with h5py.File(path, 'r') as f:
        v = np.asarray(f['mean_lr_sum'], dtype=np.float32).transpose(2, 1, 0)
        m = np.asarray(f['brainMask_merged'], dtype=np.uint8).transpose(2, 1, 0) > 0
    return v, m


def block2(a, how='mean'):
    """Downsample by 2 in every dimension (10 um -> 20 um)."""
    s = [d // 2 * 2 for d in a.shape]
    a = a[:s[0], :s[1], :s[2]]
    b = a.reshape(s[0] // 2, 2, s[1] // 2, 2, s[2] // 2, 2)
    if how == 'mean':
        return b.mean(axis=(1, 3, 5), dtype=np.float32)
    return b[:, 0, :, 0, :, 0]          # nearest: top-left voxel of each block


def main():
    t0 = time.time()

    # ---------------------------------------------------------------- P20
    p20, p20_mask = load_mat(P20_MAT)                          # (994, 800, 570)
    print(f'P20 mean LR-sum {p20.shape}, brain voxels {p20_mask.sum():,}', flush=True)
    # Only voxels where EVERY mouse has real tissue in both hemispheres. Where a
    # brain's sections have run out, P6bis's background detector returns an
    # all-tissue mask for the empty plane, the zeros are normalised to
    # -intercept/slope, and P8's abs() turns that into a large positive
    # constant in the cohort mean -- a flat slab. coverage_*.npz (from the raw
    # 4D volumes and the per-mouse masks) says where that cannot have happened.
    cov = np.load(os.path.join(OUT, 'coverage_P20.npz'))['full20']           # (497, 400, 285)
    p20_20 = block2(np.where(p20_mask, p20, 0.0))               # (497, 400, 285)
    p20_m20 = (block2(p20_mask.astype(np.float32)) > 0.5) & cov
    p20_20[~p20_m20] = 0.0
    print(f'  fully covered by all 3 mice: {cov.sum():,} of {(block2(p20_mask.astype(np.float32)) > 0.5).sum():,} 20 um brain voxels', flush=True)

    # back onto the full DeMBA 20 um grid (705 x 400 x 570), both hemispheres
    full = np.zeros((705, 400, 570), np.float32)
    fmask = np.zeros((705, 400, 570), np.uint8)
    a0, a1 = DEMBA_APLIMS[0] - 1, DEMBA_APLIMS[1]                # 0-based slice
    assert a1 - a0 == p20_20.shape[0], (a1 - a0, p20_20.shape)
    full[a0:a1, :, :285] = p20_20
    full[a0:a1, :, 285:] = p20_20[:, :, ::-1]
    fmask[a0:a1, :, :285] = p20_m20
    fmask[a0:a1, :, 285:] = p20_m20[:, :, ::-1]

    # DeMBA -> CCF, P20 -> P56. demba_dev_mouse wants (ML, DV, AP).
    print('transforming P20 LR-sum into the adult CCF ...', flush=True)
    vol = Volume(values=np.transpose(full, (2, 1, 0)), space='demba_dev_mouse',
                 voxel_size_micron=20, age_PND=20, segmentation_file=False)
    vol.transform(target_age=56, target_space='allen_mouse')
    p20_ccf = np.asarray(vol.values, dtype=np.float32)          # (AP, DV, ML) 20 um
    vm = Volume(values=np.transpose(fmask, (2, 1, 0)), space='demba_dev_mouse',
                voxel_size_micron=20, age_PND=20, segmentation_file=True)
    vm.transform(target_age=56, target_space='allen_mouse')
    p20_ccf_mask = np.asarray(vm.values) > 0
    print(f'  -> {p20_ccf.shape} in allen_mouse, {time.time() - t0:.0f} s', flush=True)
    nib.save(nib.Nifti1Image(p20_ccf, np.diag([0.02, 0.02, 0.02, 1])),
             os.path.join(OUT, 'P20_lrsum_in_ccf20.nii.gz'))

    # crop to the adult analysis grid and keep the left half
    c0 = CCF_APLIMS[0] // 2                                       # 90: 20 um voxel fully inside the crop
    c1 = c0 + (CCF_APLIMS[1] - CCF_APLIMS[0] + 1) // 2            # 540
    p20_h = p20_ccf[c0:c1, :, :285]
    p20_hm = p20_ccf_mask[c0:c1, :, :285]

    # -------------------------------------------------------------- adult
    adu, adu_mask = load_mat(ADU_MAT)                            # (900, 800, 570)
    print(f'adult mean LR-sum {adu.shape}, brain voxels {adu_mask.sum():,}', flush=True)
    adu_h = block2(np.where(adu_mask, adu, 0.0))                 # (450, 400, 285)
    adu_hm = block2(adu_mask.astype(np.float32)) > 0.5
    # same rule for the adults: all five naive and all five rws brains present
    cov_a = (np.load(os.path.join(OUT, 'coverage_naive.npz'))['full20'] &
             np.load(os.path.join(OUT, 'coverage_rws.npz'))['full20'])
    print(f'  adults fully covered by all 10 mice: {cov_a.sum():,} of {adu_hm.sum():,} 20 um brain voxels', flush=True)
    adu_hm &= cov_a
    adu_h[~adu_hm] = 0.0
    assert adu_h.shape == p20_h.shape, (adu_h.shape, p20_h.shape)

    # ------------------------------------------------------ atlas labels
    ann = np.asarray(nib.load(CCF_ANNOT).dataobj)                 # (1320, 800, 1140) 10 um
    ann = ann[CCF_APLIMS[0] - 1:CCF_APLIMS[1]]                    # 900 planes
    ann20 = block2(ann, 'nearest')[:, :, :285]                    # (450, 400, 285)
    del ann
    inside = ann20 > 0

    # ------------------------------------------- relative levels + ratio
    r0, r1 = ANALYSIS_10UM[0] // 2, ANALYSIS_10UM[1] // 2
    arange = np.zeros_like(inside)
    arange[r0:r1] = True
    both = inside & adu_hm & p20_hm & arange

    # Scale each cohort by the MEDIAN of its own isocortex, not the brain mean.
    # P6bis already anchored every mouse on the cohort's cortex, so this is the
    # scaling consistent with the pipeline; and a median is not dragged around
    # by the hippocampus, which in P20 is several times the rest of the brain
    # and would otherwise depress every cortical number. The question this
    # answers is "which areas hold a larger share of the CORTICAL signal in
    # pups than in adults" -- the relative form of the critical-period
    # prediction. It cannot say whether pups have more label in absolute terms;
    # nothing in these volumes can, after P6bis.
    iso_ids = set()
    with open(CSV_MAP, newline='', encoding='utf-8') as fh:
        for row in csv.DictReader(fh):
            if row['parcellation_term_set_name'] == 'division' and row['parcellation_term_acronym'] == 'Isocortex':
                iso_ids.add(int(row['parcellation_index']))
    iso = both & np.isin(ann20, list(iso_ids))
    adu_rel = adu_h / np.median(adu_h[iso])
    p20_rel = p20_h / np.median(p20_h[iso])
    print(f'isocortex median used for scaling -- adult {np.median(adu_h[iso]):.3f}, P20 {np.median(p20_h[iso]):.3f} '
          f'(cohort units, not comparable) over {iso.sum():,} cortical and {both.sum():,} shared voxels', flush=True)

    sm = lambda v: gaussian_filter(v, SMOOTH_SIGMA_20UM)
    eps = 1e-3
    log2r = np.log2((sm(p20_rel) + eps) / (sm(adu_rel) + eps))
    log2r[~both] = np.nan

    np.savez_compressed(os.path.join(OUT, 'volumes_ccf20.npz'), adu_rel=adu_rel.astype(np.float32),
                        p20_rel=p20_rel.astype(np.float32), log2r=log2r.astype(np.float32),
                        both=both, annot20=ann20.astype(np.uint16))

    # ----------------------------------------------------- per-structure
    names, acro, divi = {}, {}, {}
    with open(CSV_MAP, newline='', encoding='utf-8') as fh:
        for row in csv.DictReader(fh):
            idx = int(row['parcellation_index'])
            if row['parcellation_term_set_name'] == 'structure':
                names[idx] = row['parcellation_term_name']; acro[idx] = row['parcellation_term_acronym']
            elif row['parcellation_term_set_name'] == 'division':
                divi[idx] = row['parcellation_term_acronym']
    labs = ann20[both]
    sa, sp, cnt = defaultdict(float), defaultdict(float), defaultdict(int)
    for lab, a, p in zip(labs, adu_rel[both], p20_rel[both]):
        sa[lab] += a; sp[lab] += p; cnt[lab] += 1
    rows = []
    for lab in cnt:
        if cnt[lab] < 200:            # at 20 um: 1.6 nl -- skip slivers
            continue
        ma, mp = sa[lab] / cnt[lab], sp[lab] / cnt[lab]
        rows.append((names.get(lab, f'id{lab}'), acro.get(lab, ''), divi.get(lab, ''),
                     ma, mp, np.log2((mp + eps) / (ma + eps)), cnt[lab]))
    rows.sort(key=lambda r: -r[5])
    with open(os.path.join(OUT, 'region_log2ratio.csv'), 'w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(['structure', 'acronym', 'division', 'adult_rel', 'P20_rel', 'log2_P20_over_adult', 'voxels_20um'])
        for r in rows:
            w.writerow([r[0], r[1], r[2], f'{r[3]:.4f}', f'{r[4]:.4f}', f'{r[5]:.3f}', r[6]])
    print(f'\n{len(rows)} structures with >=200 shared voxels. Relative level = fraction of own brain mean.')
    print(f'{"":36s} {"adult":>7s} {"P20":>7s} {"log2":>6s}')
    print('MORE nanobody share in P20:')
    for r in rows[:15]:
        print(f'  {r[1]:8s} {r[0][:26]:26s} {r[3]:7.2f} {r[4]:7.2f} {r[5]:+6.2f}   [{r[2]}]')
    print('LESS nanobody share in P20:')
    for r in rows[-15:]:
        print(f'  {r[1]:8s} {r[0][:26]:26s} {r[3]:7.2f} {r[4]:7.2f} {r[5]:+6.2f}   [{r[2]}]')
    # The question actually asked: within the cortex, do visual and
    # somatosensory areas hold a larger share in pups?
    print('\nISOCORTEX, as a share of each cohort\'s own cortex median (1.00 = a typical cortical voxel):')
    cortex = [r for r in rows if r[2] == 'Isocortex']
    for r in sorted(cortex, key=lambda r: -r[5]):
        tag = '  <-- visual' if r[1].startswith('VIS') else ('  <-- somatosensory' if r[1].startswith('SS') else '')
        print(f'  {r[1]:9s} {r[0][:34]:34s} {r[3]:6.2f} {r[4]:6.2f} {r[5]:+6.2f}{tag}')

    # --------------------------------------------------------------- figure
    planes = [150, 250, 320]                                       # 20 um, adult crop: ~10 um 300/500/640
    vmax = np.nanpercentile(adu_rel[both], 99)
    fig, axes = plt.subplots(len(planes), 3, figsize=(15, 4.6 * len(planes)))
    for i, z in enumerate(planes):
        for j, (im, title, cmap, lim) in enumerate((
                (adu_rel[z], 'adult (naive + rws, n = 10)', 'hot', (0, vmax)),
                (p20_rel[z], 'P20 in CCF (n = 3)', 'hot', (0, vmax)),
                (log2r[z], 'log2( P20 / adult )', 'RdBu_r', (-2, 2)))):
            ax = axes[i, j]
            show = np.where(inside[z], im, np.nan) if j < 2 else im
            h = ax.imshow(np.rot90(np.where(np.isfinite(show), show, np.nan), -1)[:, ::-1] if False else show.T,
                          cmap=cmap, vmin=lim[0], vmax=lim[1], origin='lower', interpolation='nearest')
            ax.contour(ann20[z].T > 0, levels=[0.5], colors='w', linewidths=0.3)
            ax.set_title(f'{title}   plane {z * 2 + CCF_APLIMS[0]} (10 um)', fontsize=10)
            ax.set_xticks([]); ax.set_yticks([])
            plt.colorbar(h, ax=ax, fraction=0.03, pad=0.01)
    fig.suptitle('Nano LR-sum, each cohort scaled to its own brain mean, P20 carried into the adult CCF '
                 '(DeMBA -> CCF transform)', fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'slices_adult_vs_P20.png'), dpi=110)
    print(f'\nwrote {OUT}  ({time.time() - t0:.0f} s)')


if __name__ == '__main__':
    main()
