"""
v2, step 3: P20 vs adult voxelwise, in the adult CCF.

Takes the cohort volumes from v2_cohort.py (20 um, own atlas), averages the
two hemispheres, carries the P20 volumes into CCF with the DeMBA -> CCF
transform (brainglobe_ccf_translator, the same call that was validated on
the annotation), and compares:

  ratio   log2( P20 nano/auto  /  adult nano/auto )     absolute-ish
  cref    log2( P20 cortex-relative / adult cortex-relative )   distribution

Each voxel needs at least MIN_N_P20 pups and MIN_N_ADULT adults with tissue
there; the N maps travel through the transform too (nearest), so the
threshold is applied in CCF. Grey in the figure = below that.

Outputs, in data/comparisons_v2/young_P20_vs_adult/:
  volumes_ccf20.npz         adult_*, p20_*, log2_*, n_*, annot20 (AP, DV, ML-half)
  slices_*.png              one figure per mode, dorsal-up, midline right
  region_table.csv          per structure, both modes, from the CCF volumes
  cortex_table.txt          the cortical areas, printed

  D:\\sep_histology\\code\\tools\\venv_atlas\\Scripts\\python.exe v2_compare.py
"""

import csv
import os
import time
from collections import defaultdict

import numpy as np
import nibabel as nib
from scipy.ndimage import gaussian_filter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DATA = r'D:\sep_histology\data'
V2 = os.path.join(DATA, 'comparisons_v2')
OUT = os.path.join(V2, 'young_P20_vs_adult')
CSV_MAP = os.path.join(DATA, 'atlas', 'parcellation_to_parcellation_term_membership.csv')
CCF_AP0, CCF_AP1 = 180, 1079         # registered adult crop, 10 um planes
DEMBA_AP0 = 62                       # registered P20 crop starts at 20 um plane 62 of the 705-plane DeMBA volume
DEMBA_SHAPE = (705, 400, 570)        # DeMBA P20 at 20 um, (AP, DV, ML)
MIN_N_P20, MIN_N_ADULT = 2, 5
SMOOTH = 1.0                         # voxels at 20 um, on the log2 map only
MODES = ('ratio', 'cref')


def fold(v):
    """Average the two hemispheres of an (AP, DV, ML) volume -> (AP, DV, ML/2), NaN-aware."""
    ml = v.shape[2]; h = ml // 2
    left, right = v[:, :, :h], v[:, :, ml - h:][:, :, ::-1]
    return np.nanmean(np.stack([left, right]), axis=0)


def fold_n(n):
    ml = n.shape[2]; h = ml // 2
    return np.maximum(n[:, :, :h], n[:, :, ml - h:][:, :, ::-1])


def to_ccf(vol_demba_full, is_mask=False):
    """(AP, DV, ML) DeMBA P20 20 um volume -> (AP, DV, ML) CCF 20 um volume, via the translator."""
    from brainglobe_ccf_translator import Volume
    # demba_dev_mouse wants (ML, DV, AP); transform() works in place and the
    # result comes back already as (AP, DV, ML) -- the same call that was
    # validated on the annotation in compare_young_vs_adult_lrsum.py
    v = Volume(values=np.ascontiguousarray(np.transpose(vol_demba_full, (2, 1, 0))),
               space='demba_dev_mouse', voxel_size_micron=20, age_PND=20, segmentation_file=is_mask)
    v.transform(target_age=56, target_space='allen_mouse')
    return np.asarray(v.values, dtype=np.float32)


def unfold_into_demba(half):
    """(AP_crop, DV, ML/2) -> full (705, 400, 570) DeMBA grid, mirrored across the midline."""
    full = np.full(DEMBA_SHAPE, np.nan, np.float32)
    ap = half.shape[0]; h = half.shape[2]
    full[DEMBA_AP0:DEMBA_AP0 + ap, :, :h] = half
    full[DEMBA_AP0:DEMBA_AP0 + ap, :, 570 - h:] = half[:, :, ::-1]
    return full


def draw_figures(z, out):
    """Slices figure per mode from the saved volumes (also used by v2_replot.py).

    Colour conventions, so that nothing is ambiguous: no data (outside the
    atlas, or too few mice) is flat GREY; the intensity map is 'hot' cut off
    before its white end, so saturation is bright yellow and never white; the
    range is the 99th percentile of the ADULT ISOCORTEX, the structure the
    question is about, so the hippocampus saturates by design.
    """
    from matplotlib.colors import LinearSegmentedColormap
    ann_h, both = z['annot20'], z['both']
    inside = ann_h > 0
    hot = plt.get_cmap('hot')
    hot_cut = LinearSegmentedColormap.from_list('hot_cut', hot(np.linspace(0, 0.82, 256)))
    grey = '#bfbfbf'
    iso_ids = set()
    with open(CSV_MAP, newline='', encoding='utf-8') as fh:
        for row in csv.DictReader(fh):
            if row['parcellation_term_set_name'] == 'division' and row['parcellation_term_acronym'] == 'Isocortex':
                iso_ids.add(int(row['parcellation_index']))
    iso = both & np.isin(ann_h, list(iso_ids))
    cov = both.reshape(both.shape[0], -1).sum(1)
    ok = np.nonzero(cov > 0.5 * cov.max())[0]
    planes = [int(v) for v in np.linspace(ok[0], ok[-1], 6).round()]
    what = {'ratio': 'nano / autofluorescence, both background-subtracted',
            'cref': "background-subtracted nano relative to each mouse's isocortex mean"}
    for m in MODES:
        adult_v, p20_v, log2_v = z[f'adult_{m}'], z[f'p20_{m}'], z[f'log2_{m}']
        vmax = float(np.nanpercentile(adult_v[iso], 99))
        lim2 = float(np.nanpercentile(np.abs(log2_v[both]), 98))
        fig, axes = plt.subplots(len(planes), 3, figsize=(13.5, 3.9 * len(planes)))
        for i, zc in enumerate(planes):
            shown = both[zc]
            panels = ((np.where(inside[zc] & np.isfinite(adult_v[zc]), adult_v[zc], np.nan), f'adult (n = 10)  {m}', hot_cut, (0, vmax)),
                      (np.where(shown, p20_v[zc], np.nan), f'P20 in CCF (n = 3)  {m}', hot_cut, (0, vmax)),
                      (np.where(shown, log2_v[zc], np.nan), 'log2( P20 / adult )', 'RdBu_r', (-lim2, lim2)))
            for j, (im, title, cmap, lim) in enumerate(panels):
                ax = axes[i, j]
                bg = np.zeros(inside[zc].shape + (4,)); bg[inside[zc]] = matplotlib.colors.to_rgba(grey)
                ax.imshow(bg, origin='upper', interpolation='nearest', aspect='equal')
                h = ax.imshow(im, cmap=cmap, vmin=lim[0], vmax=lim[1], origin='upper', interpolation='nearest', aspect='equal')
                ax.set_title(f'{title}   plane {zc * 2 + CCF_AP0} / 10 um', fontsize=9.5)
                ax.set_xticks([]); ax.set_yticks([])
                plt.colorbar(h, ax=ax, fraction=0.035, pad=0.01, extend='max' if j < 2 else 'both')
        fig.suptitle(f'{what[m]}.  Hemispheres averaged, P20 carried DeMBA -> CCF.\n'
                     f'Grey = no data (fewer than {MIN_N_P20} pups or {MIN_N_ADULT} adults with tissue).  '
                     f'Colour range 0 to {vmax:.2f} = 99th percentile of adult isocortex; above that is bright yellow, never white.',
                     fontsize=10)
        fig.tight_layout(rect=(0, 0, 1, 0.975))
        fig.savefig(os.path.join(out, f'slices_{m}.png'), dpi=105)
        plt.close(fig)


def main():
    t0 = time.time()
    os.makedirs(OUT, exist_ok=True)
    # 20 um CCF planes 90..539 = 10 um planes 180..1079, the adult registered crop
    ann = np.asarray(nib.load(os.path.join(DATA, 'atlas', 'annotation_10.nii.gz')).dataobj)[CCF_AP0:CCF_AP1 + 1:2, ::2, ::2]
    ann_h = ann[:, :, :285]                                              # left half, (450, 400, 285)
    inside = ann_h > 0

    adult = {m: fold(np.load(os.path.join(V2, 'adult', f'{m}_mean.npy'))) for m in MODES}
    adult_n = fold_n(np.load(os.path.join(V2, 'adult', 'cref_n.npy')))
    p20 = {m: fold(np.load(os.path.join(V2, 'young_P20', f'{m}_mean.npy'))) for m in MODES}
    p20_n = fold_n(np.load(os.path.join(V2, 'young_P20', 'cref_n.npy')))
    print(f'loaded, adult {adult["cref"].shape}, P20 {p20["cref"].shape}   {time.time() - t0:.0f} s', flush=True)

    # P20 -> CCF. The translator wants finite values; carry a validity mask and
    # the n map along, and re-mask afterwards.
    p20_ccf, ok_ccf = {}, None
    valid = np.isfinite(p20['cref']) & (p20_n >= MIN_N_P20)
    ok_ccf = to_ccf(np.nan_to_num(unfold_into_demba(valid.astype(np.float32)), nan=0.0)) > 0.5
    n_ccf = np.rint(to_ccf(np.nan_to_num(unfold_into_demba(p20_n.astype(np.float32)), nan=0.0), is_mask=True)).astype(np.int16)
    for m in MODES:
        p20_ccf[m] = to_ccf(np.nan_to_num(unfold_into_demba(np.where(valid, p20[m], 0).astype(np.float32)), nan=0.0))
        print(f'{m}: P20 carried into CCF   {time.time() - t0:.0f} s', flush=True)
    # crop the CCF result to the adult registered grid, left half
    sl = slice(CCF_AP0 // 2, CCF_AP0 // 2 + ann.shape[0])
    ok_h = ok_ccf[sl, :, :285]; n_h = n_ccf[sl, :, :285]
    p20_h = {m: np.where(ok_h, p20_ccf[m][sl, :, :285], np.nan) for m in MODES}

    both = inside & ok_h & (adult_n >= MIN_N_ADULT) & np.isfinite(adult['cref'])
    print(f'voxels compared: {both.sum():,} of {inside.sum():,} inside the atlas '
          f'(adult n>={MIN_N_ADULT}: {(inside & (adult_n >= MIN_N_ADULT)).sum():,}; P20 n>={MIN_N_P20}: {(inside & ok_h).sum():,})', flush=True)

    eps = {'ratio': 0.02, 'cref': 0.02}
    log2 = {}
    for m in MODES:
        a = np.where(both, adult[m], np.nan); p = np.where(both, p20_h[m], np.nan)
        r = np.log2(np.maximum(p, eps[m]) / np.maximum(a, eps[m]))
        # smooth the map, NaN-aware
        w = both.astype(np.float32)
        num = gaussian_filter(np.nan_to_num(r) * w, SMOOTH); den = gaussian_filter(w, SMOOTH)
        log2[m] = np.where(both, num / np.maximum(den, 1e-3), np.nan).astype(np.float32)
    np.savez_compressed(os.path.join(OUT, 'volumes_ccf20.npz'), annot20=ann_h, both=both, adult_n=adult_n, p20_n=n_h,
                        **{f'adult_{m}': adult[m].astype(np.float32) for m in MODES},
                        **{f'p20_{m}': p20_h[m].astype(np.float32) for m in MODES},
                        **{f'log2_{m}': log2[m] for m in MODES})

    # ------------------------------------------------ per-structure table
    names, acro, divi = {}, {}, {}
    with open(CSV_MAP, newline='', encoding='utf-8') as fh:
        for row in csv.DictReader(fh):
            idx = int(row['parcellation_index'])
            if row['parcellation_term_set_name'] == 'structure':
                names[idx] = row['parcellation_term_name']; acro[idx] = row['parcellation_term_acronym']
            elif row['parcellation_term_set_name'] == 'division':
                divi[idx] = row['parcellation_term_acronym']
    labs = ann_h[both]
    acc = {m: (defaultdict(float), defaultdict(float)) for m in MODES}
    cnt = defaultdict(int); meta = {}
    for m in MODES:
        for lab, a, p in zip(labs, adult[m][both], p20_h[m][both]):
            key = names.get(int(lab), f'id{lab}')
            acc[m][0][key] += a; acc[m][1][key] += p
            if m == 'cref':
                cnt[key] += 1; meta[key] = (acro.get(int(lab), ''), divi.get(int(lab), ''))
    rows = []
    for key in cnt:
        if cnt[key] < 100:
            continue
        r = {'structure': key, 'acronym': meta[key][0], 'division': meta[key][1], 'voxels_20um': cnt[key]}
        for m in MODES:
            a, p = acc[m][0][key] / cnt[key], acc[m][1][key] / cnt[key]
            r[f'adult_{m}'] = a; r[f'P20_{m}'] = p
            r[f'log2_{m}'] = np.log2(max(p, eps[m]) / max(a, eps[m]))
        rows.append(r)
    rows.sort(key=lambda r: (r['division'], r['acronym']))
    with open(os.path.join(OUT, 'region_table.csv'), 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader()
        for r in rows:
            w.writerow({k: (f'{v:.4f}' if isinstance(v, float) else v) for k, v in r.items()})
    cortex = sorted([r for r in rows if r['division'] == 'Isocortex'], key=lambda r: -r['log2_cref'])
    lines = [f'{"area":9s} {"structure":34s} {"adult r":>8s} {"P20 r":>7s} {"log2 ratio":>10s} | {"adult c":>8s} {"P20 c":>7s} {"log2 cref":>9s}']
    for r in cortex:
        tag = '  <-- visual' if r['acronym'].startswith('VIS') else ('  <-- somatosensory' if r['acronym'].startswith('SS') else '')
        lines.append(f'{r["acronym"]:9s} {r["structure"][:34]:34s} {r["adult_ratio"]:8.2f} {r["P20_ratio"]:7.2f} {r["log2_ratio"]:+10.2f} | '
                     f'{r["adult_cref"]:8.2f} {r["P20_cref"]:7.2f} {r["log2_cref"]:+9.2f}{tag}')
    with open(os.path.join(OUT, 'cortex_table.txt'), 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    print('\n'.join(lines))

    draw_figures(dict(annot20=ann_h, both=both, adult_n=adult_n, p20_n=n_h,
                      **{f'adult_{m}': adult[m] for m in MODES}, **{f'p20_{m}': p20_h[m] for m in MODES},
                      **{f'log2_{m}': log2[m] for m in MODES}), OUT)
    print(f'done  {time.time() - t0:.0f} s')


if __name__ == '__main__':
    main()
