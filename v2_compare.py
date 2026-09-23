"""
v2, step 4: young vs adult in the adult CCF, from cohort volumes that are
already there.

Nothing is transformed here any more: v2_to_ccf carried each brain onto the
660 x 400 x 570 grid at 20 um and v2_cohort averaged them on it. This script
folds the hemispheres, compares, and writes the maps and tables.

  ratio   log2( young nano/auto  /  adult nano/auto )        absolute-ish
  cref    log2( young cortex-relative / adult cortex-relative )  distribution

The young group is P20 + P16 pooled (YOUNG below); young_P20 is written too, as
the sensitivity check for what the single P16 brain does to the answer. A voxel
enters the comparison only where at least MIN_N_YOUNG young and MIN_N_ADULT
adult brains have tissue; the rest is grey in the figures.

Outputs, in data/comparisons_v2/young_vs_adult/:
  volumes_ccf20.npz     adult_*, young_*, log2_*, n maps, annot20 (AP, DV, ML-half)
  slices_<mode>.png     dorsal-up, midline right, no data in grey
  region_table.csv      per structure, both modes, both young groups
  cortex_table.txt      the cortical areas, printed

  D:\\sep_histology\\code\\tools\\venv_atlas\\Scripts\\python.exe v2_compare.py
"""

import csv
import os
import time

import numpy as np
import nibabel as nib
from scipy.ndimage import gaussian_filter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from v2_per_mouse import DATA, CSV_MAP
from v2_cohort import OUT_ROOT as CCF_ROOT, COHORTS, MODES

OUT = os.path.join(DATA, 'comparisons_v2', 'young_vs_adult')
CCF_AP0 = 180                        # the adult registered crop, 10 um planes
MIN_N_YOUNG, MIN_N_ADULT = 2, 5
SMOOTH = 1.0                         # voxels at 20 um, applied to the log2 map only
from v2_cohort import MODES
YOUNG = 'young'                      # pooled P20 + P16
YOUNG_ALT = 'young_P20'              # sensitivity check


def save_figure(fig, path):
    """Save, and if the file is open in a viewer say so instead of dying.

    Windows refuses to overwrite a PNG that an image viewer holds open, and a
    run that writes several figures should not lose the rest because one of
    them was being looked at.
    """
    try:
        fig.savefig(path, dpi=105)
    except OSError:
        alt = path.replace('.png', '_new.png')
        fig.savefig(alt, dpi=105)
        print(f'  NOTE: {os.path.basename(path)} is open elsewhere; wrote {os.path.basename(alt)} instead', flush=True)


def fold(v):
    """Average the two hemispheres of an (AP, DV, ML) volume -> (AP, DV, ML/2), NaN-aware."""
    ml = v.shape[2]; h = ml // 2
    left, right = v[:, :, :h], v[:, :, ml - h:][:, :, ::-1]
    return np.nanmean(np.stack([left, right]), axis=0)


def fold_n(n):
    ml = n.shape[2]; h = ml // 2
    return np.maximum(n[:, :, :h], n[:, :, ml - h:][:, :, ::-1])


def load_cohort(cohort):
    """(modes dict folded, folded n map) for a cohort in CCF."""
    m = {k: fold(np.load(os.path.join(CCF_ROOT, cohort, f'{k}_mean.npy'))) for k in MODES}
    n = fold_n(np.load(os.path.join(CCF_ROOT, cohort, 'cref_n.npy')))
    return m, n


def draw_figures(z, out):
    """Slices figure per mode from the saved volumes (also used by v2_replot.py).

    Colour conventions: no data (outside the atlas, or too few mice) is flat
    GREY; the intensity map stops before the white end of 'hot', so saturation
    reads as bright yellow and is never confused with empty; the range is the
    99th percentile of the ADULT ISOCORTEX, the structure the question is
    about, so the hippocampus saturates by design.
    """
    from matplotlib.colors import LinearSegmentedColormap
    ann_h, both = z['annot20'], z['both']
    inside = ann_h > 0
    hot = plt.get_cmap('hot')
    hot_cut = LinearSegmentedColormap.from_list('hot_cut', hot(np.linspace(0, 0.82, 256)))
    hot_cut.set_bad((0, 0, 0, 0))
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
            'cref': "background-subtracted nano relative to each mouse's isocortex mean",
            'subref': 'background-subtracted nano relative to the subcortex, excluding HPF and STR',
            'zref': "range-matched: position within each brain's own distribution "
                    '(median 0, p90-p10 = 1), so level and dynamic range are equal across brains'}
    n_young = int(z['n_young_mice']); n_adult = int(z['n_adult_mice'])
    for m in MODES:
        adult_v, young_v, log2_v = z[f'adult_{m}'], z[f'young_{m}'], z[f'log2_{m}']
        signed = m == 'zref'
        vmax = float(np.nanpercentile(adult_v[iso], 99)) if not signed else \
            float(np.nanpercentile(np.abs(adult_v[both]), 98))
        lim2 = float(np.nanpercentile(np.abs(log2_v[both]), 98))
        fig, axes = plt.subplots(len(planes), 3, figsize=(13.5, 3.9 * len(planes)))
        for i, zc in enumerate(planes):
            shown = both[zc]
            cmap_mean = 'RdBu_r' if signed else hot_cut
            lim_mean = (-vmax, vmax) if signed else (0, vmax)
            diff_name = 'young - adult' if signed else 'log2( young / adult )'
            panels = ((np.where(inside[zc] & np.isfinite(adult_v[zc]), adult_v[zc], np.nan),
                       f'adult (n = {n_adult})  {m}', cmap_mean, lim_mean),
                      (np.where(shown, young_v[zc], np.nan), f'young (n = {n_young})  {m}', cmap_mean, lim_mean),
                      (np.where(shown, log2_v[zc], np.nan), diff_name, 'RdBu_r', (-lim2, lim2)))
            for j, (im, title, cmap, lim) in enumerate(panels):
                ax = axes[i, j]
                bg = np.zeros(inside[zc].shape + (4,)); bg[inside[zc]] = matplotlib.colors.to_rgba(grey)
                ax.imshow(bg, origin='upper', interpolation='nearest', aspect='equal')
                h = ax.imshow(im, cmap=cmap, vmin=lim[0], vmax=lim[1], origin='upper',
                              interpolation='nearest', aspect='equal')
                ax.set_title(f'{title}   plane {zc * 2 + CCF_AP0} / 10 um', fontsize=9.5)
                ax.set_xticks([]); ax.set_yticks([])
                plt.colorbar(h, ax=ax, fraction=0.035, pad=0.01, extend='max' if j < 2 else 'both')
        scale_note = (f"diverging scale, 0 = that brain's median structure, +-{vmax:.2f} = its own p10-p90 spread"
                      if signed else
                      f'colour range 0 to {vmax:.2f} = 99th percentile of adult isocortex; brighter is yellow, never white')
        fig.suptitle('\n'.join((what[m] + '.',
                                'Hemispheres averaged; every brain carried into the adult CCF individually.  '
                                f'Grey = no data (fewer than {MIN_N_YOUNG} young or {MIN_N_ADULT} adults with tissue).',
                                scale_note)), fontsize=10)
        fig.tight_layout(rect=(0, 0, 1, 0.975))
        save_figure(fig, os.path.join(out, f'slices_{m}.png'))
        plt.close(fig)


def main():
    t0 = time.time()
    os.makedirs(OUT, exist_ok=True)
    # 20 um CCF planes 90..539 hold the adult registered crop; the young volumes
    # live on the same 660-plane grid, so the annotation is simply halved here.
    ann = np.asarray(nib.load(os.path.join(DATA, 'atlas', 'annotation_10.nii.gz')).dataobj)[::2, ::2, ::2]
    ann_h = ann[:, :, :285]
    inside = ann_h > 0

    adult, adult_n = load_cohort('adult')
    young, young_n = load_cohort(YOUNG)
    young_alt, young_alt_n = load_cohort(YOUNG_ALT)
    both = inside & (young_n >= MIN_N_YOUNG) & (adult_n >= MIN_N_ADULT) & np.isfinite(young['cref']) & np.isfinite(adult['cref'])
    print(f'voxels compared: {both.sum():,} of {inside.sum():,} inside the atlas '
          f'(young n>={MIN_N_YOUNG}: {(inside & (young_n >= MIN_N_YOUNG)).sum():,}; '
          f'adult n>={MIN_N_ADULT}: {(inside & (adult_n >= MIN_N_ADULT)).sum():,})   {time.time() - t0:.0f} s', flush=True)

    eps = 0.02
    log2, log2_alt = {}, {}
    for m in MODES:
        for src, dst in ((young, log2), (young_alt, log2_alt)):
            a = np.where(both, adult[m], np.nan); p = np.where(both, src[m], np.nan)
            # zref is already a log-scale position, so the two groups are
            # compared by difference; everything else by log2 ratio
            r = (p - a) if m == 'zref' else np.log2(np.maximum(p, eps) / np.maximum(a, eps))
            w = both.astype(np.float32)
            num = gaussian_filter(np.nan_to_num(r) * w, SMOOTH); den = gaussian_filter(w, SMOOTH)
            dst[m] = np.where(both, num / np.maximum(den, 1e-3), np.nan).astype(np.float32)
    np.savez_compressed(os.path.join(OUT, 'volumes_ccf20.npz'), annot20=ann_h, both=both,
                        adult_n=adult_n, young_n=young_n, young_alt_n=young_alt_n,
                        n_young_mice=len(COHORTS[YOUNG]), n_adult_mice=len(COHORTS['adult']),
                        **{f'adult_{m}': adult[m].astype(np.float32) for m in MODES},
                        **{f'young_{m}': young[m].astype(np.float32) for m in MODES},
                        **{f'young_alt_{m}': young_alt[m].astype(np.float32) for m in MODES},
                        **{f'log2_{m}': log2[m] for m in MODES},
                        **{f'log2_alt_{m}': log2_alt[m] for m in MODES})

    # ------------------------------------------------ per-structure table
    names, acro, divi = {}, {}, {}
    with open(CSV_MAP, newline='', encoding='utf-8') as fh:
        for row in csv.DictReader(fh):
            idx = int(row['parcellation_index'])
            if row['parcellation_term_set_name'] == 'structure':
                names[idx] = row['parcellation_term_name']; acro[idx] = row['parcellation_term_acronym']
            elif row['parcellation_term_set_name'] == 'division':
                divi[idx] = row['parcellation_term_acronym']
    # Aggregate per STRUCTURE, not per parcellation_index: in this ontology the
    # layers of an area are separate indices that share one structure name, and
    # a table of areas is what anyone reads. A lookup from index to structure
    # does that grouping once, and bincount then sums 18 million voxels without
    # a Python loop. A group with no value at a voxel (the P20-only group inside
    # the pooled group's mask) is left out of its own count rather than
    # poisoning its mean.
    struct_of = np.zeros(int(ann_h.max()) + 1, np.int64)
    struct_names = []
    seen = {}
    for idx in np.unique(ann_h):
        if idx == 0:
            continue
        nm = names.get(int(idx), f'id{idx}')
        if nm not in seen:
            seen[nm] = len(struct_names)
            struct_names.append((nm, acro.get(int(idx), ''), divi.get(int(idx), '')))
        struct_of[idx] = seen[nm]
    labs = struct_of[ann_h[both].astype(np.int64)]
    n_struct = len(struct_names)
    n_vox = np.bincount(labs, minlength=n_struct)
    means = {}
    for m in MODES:
        for tag, src in (('adult', adult), ('young', young), ('young_P20', young_alt)):
            v = src[m][both]
            ok = np.isfinite(v)
            tot = np.bincount(labs[ok], weights=v[ok], minlength=n_struct)
            cnt = np.bincount(labs[ok], minlength=n_struct)
            with np.errstate(invalid='ignore', divide='ignore'):
                means[(tag, m)] = np.where(cnt > 0, tot / np.maximum(cnt, 1), np.nan)

    rows = []
    for g in np.nonzero(n_vox >= 100)[0]:
        nm, ac, dv = struct_names[g]
        r = {'structure': nm, 'acronym': ac, 'division': dv, 'voxels_20um': int(n_vox[g])}
        for m in MODES:
            for tag in ('adult', 'young', 'young_P20'):
                r[f'{tag}_{m}'] = float(means[(tag, m)][g])
            if m == 'zref':
                r[f'log2_{m}'] = r[f'young_{m}'] - r[f'adult_{m}']
                r[f'log2_{m}_P20only'] = r[f'young_P20_{m}'] - r[f'adult_{m}']
            else:
                r[f'log2_{m}'] = np.log2(max(r[f'young_{m}'], eps) / max(r[f'adult_{m}'], eps))
                r[f'log2_{m}_P20only'] = np.log2(max(r[f'young_P20_{m}'], eps) / max(r[f'adult_{m}'], eps))
        rows.append(r)
    rows.sort(key=lambda r: (r['division'], r['acronym']))
    with open(os.path.join(OUT, 'region_table.csv'), 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader()
        for r in rows:
            w.writerow({k: (f'{v:.4f}' if isinstance(v, float) else v) for k, v in r.items()})

    cortex = sorted([r for r in rows if r['division'] == 'Isocortex'], key=lambda r: -r['log2_cref'])
    lines = [f'{"area":9s} {"structure":34s} ' + ' '.join(f'{m:>10s}' for m in MODES)]
    for r in cortex:
        tag = '  <-- visual' if r['acronym'].startswith('VIS') else ('  <-- somatosensory' if r['acronym'].startswith('SS') else '')
        lines.append(f'{r["acronym"]:9s} {r["structure"][:34]:34s} '
                     + ' '.join(f'{r[f"log2_{m}"]:+10.2f}' for m in MODES) + tag)
    with open(os.path.join(OUT, 'cortex_table.txt'), 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    print('\n'.join(lines))

    draw_figures(dict(np.load(os.path.join(OUT, 'volumes_ccf20.npz'))), OUT)
    print(f'done  {time.time() - t0:.0f} s')


if __name__ == '__main__':
    main()
