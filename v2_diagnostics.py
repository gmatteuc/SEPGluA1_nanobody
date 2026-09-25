"""
Diagnostic sheets for every step of the v2 route, so the pipeline can be
audited rather than believed.

One figure per question, written to
data/comparisons_v2/processing_diagnostics/, with an index README beside them:

  01_tissue_<mouse>.png       what was counted as tissue, drawn on the brain
  02_levels_<mouse>.png       where the background and the threshold sit in the
                              actual intensity distributions
  03_coverage.png             how much of the atlas each brain covers, plane by
                              plane -- i.e. where a cohort mean rests on few mice
  04_warp_<mouse>.png         the same plane before and after DeMBA -> CCF, and
                              every region's value before against after
  05_cohort_n.png             how many brains contribute at each voxel
  06_scaling.png              the per-mouse scale factors, and what they do to
                              the cortex distributions
  07_route_agreement.png      voxelwise (warped) against region-wise (never
                              warped), structure by structure
  08_mask_vs_p6bis.png        the v2 tissue mask against the P6bis/P2bis one
                              that was debugged by eye, where the latter works
  09_denominators.png         autofluorescence and SEP as reference channels,
                              and whether a structure's answer depends on which
                              one is used

The per-brain sheets P4bis writes when it carries the SEP channel into
registered space sit beside these, in processing_diagnostics/sep_channel/.

Run it after a full pass, or with mouse names to refresh a few sheets:

  D:\\sep_histology\\code\\tools\\venv_atlas\\Scripts\\python.exe v2_diagnostics.py [mouse ...]
"""

import csv
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from v2_per_mouse import annotation_20, atlas_grid, MICE, CSV_MAP, DATA, OUT as PER_MOUSE
from v2_cohort import COHORTS, OUT_ROOT as CCF_ROOT, PER_MOUSE_CCF

OUT = os.path.join(DATA, 'comparisons_v2', 'processing_diagnostics')
YOUNG = [m for m, v in MICE.items() if v[1] != 'ccf']
ADULT = [m for m, v in MICE.items() if v[1] == 'ccf']


def save_figure(fig, path, dpi=95):
    """Save, and if the file is open in a viewer say so instead of dying.

    A sheet being looked at should not cost the rest of the run.
    """
    try:
        fig.savefig(path, dpi=dpi)
    except OSError:
        alt = path.replace('.png', '_new.png')
        fig.savefig(alt, dpi=dpi)
        print(f'  NOTE: {os.path.basename(path)} is open elsewhere; wrote {os.path.basename(alt)} instead', flush=True)


def show(ax, img, mask=None, p=99.5):
    """A plane, dorsal up and ventral down, scaled to its own tissue.

    A plane of these volumes is (DV, ML), so it is drawn as it comes: rows run
    dorsal to ventral, columns left to right. The slice figures elsewhere use
    the same convention.
    """
    im = np.asarray(img, float)
    hi = np.percentile(im[im > 0], p) if (im > 0).any() else 1.0
    ax.imshow(np.clip(im / hi, 0, 1), cmap='gray', origin='upper')
    if mask is not None:
        ax.contour(np.asarray(mask, float), levels=[0.5], colors='#e74c3c', linewidths=0.9)
    ax.axis('off')


def sheet_tissue(mouse, ann, z):
    """01 -- the tissue mask on four planes, with the atlas outline for scale."""
    sig, tissue = z['sig'].astype(np.float32), z['tissue']
    brain = ann > 0
    cov = np.array([tissue[k][brain[k]].mean() if brain[k].any() else 0 for k in range(tissue.shape[0])])
    planes = np.linspace(*np.percentile(np.flatnonzero(cov > 0.2), [5, 95]), 4).round().astype(int)
    fig, axes = plt.subplots(1, 4, figsize=(19, 5.2))
    for ax, k in zip(axes, planes):
        show(ax, sig[k], tissue[k])
        ax.contour((ann[k] > 0).astype(float), levels=[0.5], colors='#3498db', linewidths=0.6)
        ax.set_title(f'plane {k}   tissue {100 * tissue[k][brain[k]].mean():.0f}% of atlas brain', fontsize=10)
    fig.suptitle(f'{mouse}: red = tissue mask (auto channel above background + 4 MAD, and a section reached here), '
                 f'blue = atlas brain.  Background subtracted: nano {float(z["bg_nano"]):.0f}, auto {float(z["bg_auto"]):.0f} counts',
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    save_figure(fig, os.path.join(OUT, f'01_tissue_{mouse}.png'), dpi=95); plt.close(fig)


def sheet_levels(mouse, ann, z):
    """02 -- the intensity distributions the mask and the background rest on."""
    sig, auto, tissue = z['sig'].astype(np.float32), z['auto'].astype(np.float32), z['tissue']
    brain = ann > 0
    bg_n, bg_a, mad = float(z['bg_nano']), float(z['bg_auto']), float(z['mad_auto'])
    off = (~brain) & (auto > -bg_a)                      # imaged, outside the atlas brain
    fig, axes = plt.subplots(1, 3, figsize=(17, 4.6))
    for ax, (arr, name, bg) in zip(axes[:2], ((auto, 'auto', bg_a), (sig, 'nano', bg_n))):
        ax.hist(arr[off][::17] + bg, bins=200, range=(0, 4 * bg), color='#95a5a6', label='off tissue', density=True)
        ax.hist(arr[tissue][::37] + bg, bins=200, range=(0, 4 * bg), color='#c0392b', alpha=0.6, label='tissue', density=True)
        ax.axvline(bg, color='k', lw=1.2, label=f'background {bg:.0f}')
        if name == 'auto':
            ax.axvline(bg + 4 * mad, color='#2980b9', lw=1.2, ls='--', label=f'mask threshold {bg + 4 * mad:.0f}')
        ax.set_xlabel(f'{name} channel, raw counts'); ax.set_ylabel('density'); ax.legend(fontsize=8)
        ax.set_title(f'{name}: tissue vs off tissue', fontsize=10)
    ax = axes[2]
    iso = tissue & np.isin(ann, ISO)
    ax.hist((sig[iso] / float(z['cortex_mean']))[::37], bins=200, range=(0, 3), color='#c0392b', density=True)
    ax.axvline(1, color='k', lw=1.2)
    ax.set_xlabel('sig / isocortex mean'); ax.set_ylabel('density')
    ax.set_title(f'cortex after scaling (mean {float(z["cortex_mean"]):.0f} counts -> 1.0)', fontsize=10)
    fig.suptitle(f'{mouse}: what the background subtraction and the mask threshold actually separate', fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    save_figure(fig, os.path.join(OUT, f'02_levels_{mouse}.png'), dpi=95); plt.close(fig)


def line_colours(n):
    """n distinguishable line colours with no green in them (plasma, trimmed)."""
    return plt.get_cmap('plasma')(np.linspace(0.0, 0.88, max(n, 2)))


def sheet_coverage():
    """03 -- per brain, the fraction of the atlas brain with tissue, plane by plane."""
    fig, axes = plt.subplots(2, 1, figsize=(14, 8))
    for ax, group, label in ((axes[0], YOUNG, 'young, on each brain\'s own atlas'),
                             (axes[1], ADULT, 'adults, on the CCF')):
        cols = line_colours(len(group))
        for ci, mouse in enumerate(group):
            ann = ANN[MICE[mouse][1]]; brain = ann > 0
            t = np.load(os.path.join(PER_MOUSE, mouse + '.npz'))['tissue']
            cov = np.array([t[k][brain[k]].mean() if brain[k].any() else np.nan for k in range(t.shape[0])])
            lo, hi = atlas_grid(MICE[mouse][1])[1]
            ax.plot(np.arange(len(cov)), 100 * cov, lw=1.3, color=cols[ci],
                    label=f'{mouse} ({np.nansum(cov > 0.2):.0f} planes)')
        ax.set_xlabel('atlas plane within the crop'); ax.set_ylabel('% of atlas brain with tissue')
        ax.set_title(label, fontsize=10); ax.grid(lw=0.3, alpha=0.6); ax.legend(fontsize=7, ncol=2)
    fig.suptitle('Coverage: where a cohort mean rests on every brain, and where it rests on one or two', fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    save_figure(fig, os.path.join(OUT, '03_coverage.png'), dpi=95); plt.close(fig)


def sheet_warp(mouse):
    """04 -- one brain before and after the transform, as a picture and as numbers."""
    ann_n = ANN[MICE[mouse][1]]
    zn = np.load(os.path.join(PER_MOUSE, mouse + '.npz'))
    zc = np.load(os.path.join(PER_MOUSE_CCF, mouse + '.npz'))
    sig_n, t_n = zn['sig'].astype(np.float32), zn['tissue']
    sig_c, t_c = zc['sig'].astype(np.float32), zc['tissue']
    # the same anatomy: match by the fraction of the covered range
    cov_n = np.flatnonzero([t_n[k].sum() > 2000 for k in range(t_n.shape[0])])
    cov_c = np.flatnonzero([t_c[k].sum() > 2000 for k in range(t_c.shape[0])])
    fig = plt.figure(figsize=(17, 5.0))
    for i, f in enumerate((0.3, 0.6)):
        kn = int(cov_n[0] + f * (cov_n[-1] - cov_n[0])); kc = int(cov_c[0] + f * (cov_c[-1] - cov_c[0]))
        ax = fig.add_subplot(1, 3, i + 1)
        # both planes are (DV, ML): stacking along ML puts them side by side
        show(ax, np.concatenate([sig_n[kn], sig_c[kc]], axis=1),
             np.concatenate([t_n[kn], t_c[kc]], axis=1))
        ax.set_title(f'{int(100 * f)}% through the stack: own atlas (left) | in CCF (right)', fontsize=10)
    ax = fig.add_subplot(1, 3, 3)
    ccf_full = np.zeros(sig_c.shape, ANN['ccf'].dtype); ccf_full[90:540] = ANN['ccf']
    a, b = [], []
    for idx in np.unique(ANN[MICE[mouse][1]])[1:][::3]:
        m1 = t_n & (ann_n == idx); m2 = t_c & (ccf_full == idx)
        if m1.sum() > 300 and m2.sum() > 300:
            a.append(sig_n[m1].mean()); b.append(sig_c[m2].mean())
    a, b = np.array(a), np.array(b)
    ax.loglog(a, b, 'o', ms=3, color='#c0392b', alpha=0.6)
    lim = [min(a.min(), b.min()) * 0.9, max(a.max(), b.max()) * 1.1]
    ax.plot(lim, lim, 'k-', lw=0.8)
    ax.set_xlabel('region mean, own atlas'); ax.set_ylabel('region mean, in CCF')
    ax.set_title(f'{len(a)} regions, median |log2 change| {np.median(np.abs(np.log2(b / a))):.3f}', fontsize=10)
    fig.suptitle(f'{mouse}: what the DeMBA -> CCF transform does (adults are placed, not warped)', fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    save_figure(fig, os.path.join(OUT, f'04_warp_{mouse}.png'), dpi=95); plt.close(fig)


def sheet_cohort_n():
    """05 -- how many brains support each voxel, per cohort."""
    cohorts = ['young', 'young_P20', 'adult']
    planes = [150, 250, 350, 450]
    fig, axes = plt.subplots(len(cohorts), len(planes), figsize=(4.2 * len(planes), 3.6 * len(cohorts)))
    for r, cohort in enumerate(cohorts):
        n = np.load(os.path.join(CCF_ROOT, cohort, 'cref_n.npy'))
        for c, k in enumerate(planes):
            ax = axes[r, c]
            h = ax.imshow(n[k], cmap='magma', vmin=0, vmax=len(COHORTS[cohort]),
                          origin='upper', interpolation='nearest')
            ax.set_title(f'{cohort}  CCF plane {2 * k} / 10 um', fontsize=9); ax.axis('off')
            plt.colorbar(h, ax=ax, fraction=0.03, pad=0.01)
    fig.suptitle('How many brains contribute at each voxel (the n map the comparison thresholds)', fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    save_figure(fig, os.path.join(OUT, '05_cohort_n.png'), dpi=95); plt.close(fig)


def sheet_scaling():
    """06 -- the per-mouse numbers every normalisation rests on."""
    rows = []
    for mouse in list(MICE):
        z = np.load(os.path.join(PER_MOUSE, mouse + '.npz'))
        ann = ANN[MICE[mouse][1]]
        iso = z['tissue'] & np.isin(ann, ISO)
        rows.append((mouse, MICE[mouse][0], float(z['bg_nano']), float(z['bg_auto']),
                     float(z['cortex_mean']), float(z['auto'].astype(np.float32)[iso].mean())))
    young = [r for r in rows if r[1].startswith('young')]
    adult = [r for r in rows if not r[1].startswith('young')]
    order = young + adult
    fig, axes = plt.subplots(1, 3, figsize=(17, 4.8))
    for ax, (j, name) in zip(axes, ((2, 'off-tissue background, nano'), (4, 'isocortex mean, nano (bg-subtracted)'))):
        for i, r in enumerate(order):
            ax.plot(i, r[j], 'o', color='#c0392b' if r[1].startswith('young') else '#555555')
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels([r[0].split('_')[0] for r in order], rotation=70, fontsize=7)
        ax.axvline(len(young) - 0.5, color='k', lw=0.6, ls=':')
        ax.set_ylabel('raw counts'); ax.set_title(name, fontsize=10); ax.grid(lw=0.3, alpha=0.6)
    ax = axes[2]
    for grp, col, lbl in ((young, '#c0392b', 'young'), (adult, '#555555', 'adult')):
        ax.plot([r[4] for r in grp], [r[5] for r in grp], 'o', color=col, label=lbl)
    ax.set_xlabel('isocortex mean, nano'); ax.set_ylabel('isocortex mean, auto')
    ax.set_title('the two channels track each other across brains', fontsize=10)
    ax.legend(fontsize=8); ax.grid(lw=0.3, alpha=0.6)
    fig.suptitle('Per-brain levels: what the background subtraction removes, and what the cortex scaling divides by\n'
                 'The 4x spread among adults is why no analysis uses raw counts', fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    save_figure(fig, os.path.join(OUT, '06_scaling.png'), dpi=95); plt.close(fig)


def sheet_route_agreement():
    """07 -- the warped voxelwise route against the never-warped region route."""
    d = os.path.join(DATA, 'comparisons_v2', 'young_vs_adult')
    vox = {r['acronym']: (float(r['log2_cref']), r['division'], int(r['voxels_20um']))
           for r in csv.DictReader(open(os.path.join(d, 'region_table.csv'), encoding='utf-8'))}
    reg = {r['acronym']: float(r['diff_log2']) for r in csv.DictReader(open(os.path.join(d, 'region_stats.csv'), encoding='utf-8'))
           if r['reading'] == 'cref'}
    keys = [a for a in vox if a in reg and vox[a][2] >= 500]
    x = np.array([vox[a][0] for a in keys]); y = np.array([reg[a] for a in keys])
    iso = np.array([vox[a][1] == 'Isocortex' for a in keys])
    fig, ax = plt.subplots(figsize=(7.2, 7))
    ax.plot(x[~iso], y[~iso], 'o', ms=4, color='#95a5a6', label=f'other ({(~iso).sum()})')
    ax.plot(x[iso], y[iso], 'o', ms=5, color='#c0392b', label=f'isocortex ({iso.sum()})')
    lim = [min(x.min(), y.min()) - 0.1, max(x.max(), y.max()) + 0.1]
    ax.plot(lim, lim, 'k-', lw=0.8); ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel('log2 young/adult -- voxelwise, every brain warped into CCF')
    ax.set_ylabel('log2 young/adult -- region-wise, no warping')
    ax.set_title(f'r = {np.corrcoef(x, y)[0, 1]:.3f} over {len(keys)} structures, '
                 f'median |difference| {np.median(np.abs(x - y)):.3f} log2\n'
                 f'isocortex only: r = {np.corrcoef(x[iso], y[iso])[0, 1]:.3f}, '
                 f'median {np.median(np.abs(x[iso] - y[iso])):.3f}', fontsize=10)
    ax.legend(fontsize=9); ax.grid(lw=0.3, alpha=0.6)
    fig.tight_layout()
    save_figure(fig, os.path.join(OUT, '07_route_agreement.png'), dpi=110); plt.close(fig)


def sheet_mask_vs_p6bis():
    """08 -- the v2 mask against the one P6bis uses, on brains that have both.

    The P6bis mask was checked by eye over many sessions, so it is the right
    thing to be measured against. It is not used here for two reasons: it is
    computed on the nano channel, which is the quantity being compared and is
    four times dimmer in pups; and its centre-of-mass guard returns an
    all-background mask for any partial section. Where it does work the two
    masks agree to better than 2% of voxels.
    """
    import h5py
    cases = [('MG903_SepGluA_P20', os.path.join(DATA, 'young', 'nano_4d_normalized_bkgmask_P20.mat'), 1, (200, 260, 330)),
             ('CGF027_Gria1', os.path.join(DATA, 'naive', 'nano_4d_normalized_bkgmask.mat'), 0, (150, 230, 300))]
    fig, axes = plt.subplots(2, 3, figsize=(18, 10.5))
    for row, (mouse, bk, idx, planes) in enumerate(cases):
        ann = ANN[MICE[mouse][1]]
        z = np.load(os.path.join(PER_MOUSE, mouse + '.npz'))
        mine, sig = z['tissue'], z['sig'].astype(np.float32)
        with h5py.File(bk, 'r') as f:
            B = f['recomputed_bkg_mask_4d']
            for col, k in enumerate(planes):
                old = np.zeros((mine.shape[1], mine.shape[2]), np.float32)
                for j in range(mine.shape[2]):
                    q = sum((np.asarray(B[idx, 2 * j + dj, :, 2 * k + dk]) == 0).astype(np.float32)
                            for dj in (0, 1) for dk in (0, 1)) / 4
                    old[:, j] = q.reshape(mine.shape[1], 2).mean(axis=1)
                ax = axes[row, col]
                show(ax, sig[k], mine[k])
                # Both are (DV, ML) like the plane underneath them, so neither is
                # transposed. The .T that used to be here belonged to the days
                # when show() turned the plane on its side.
                ax.contour((old > 0.5).astype(float), levels=[0.5], colors='#3498db', linewidths=0.9, linestyles='--')
                ax.contour((ann[k] > 0).astype(float), levels=[0.5], colors='#cccccc', linewidths=0.6)
                agree = 2 * (mine[k] & (old > 0.5)).sum() / max(mine[k].sum() + (old > 0.5).sum(), 1)
                ax.set_title(f'{mouse}  plane {k}   Dice {agree:.3f}', fontsize=10)
    fig.suptitle('red = v2 mask (auto channel), blue dashed = P6bis mask (nano channel), grey = atlas brain.\n'
                 'Swapping one for the other moves every cortical result by at most 0.02 log2', fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    save_figure(fig, os.path.join(OUT, '08_mask_vs_p6bis.png'), dpi=95); plt.close(fig)


def sheet_denominators():
    """09 -- the two reference channels side by side, and whether the answer moves."""
    rows = []
    for mouse in list(MICE):
        z = np.load(os.path.join(PER_MOUSE, mouse + '.npz'))
        if 'sep' not in z.files:
            continue
        ann = ANN[MICE[mouse][1]]
        iso = z['tissue'] & np.isin(ann, ISO)
        rows.append((mouse, MICE[mouse][0],
                     float(z['cortex_mean']),
                     float(z['auto'].astype(np.float32)[iso].mean()),
                     float(z['sep'].astype(np.float32)[iso].mean())))
    if not rows:
        print('09 skipped: no brain carries a SEP channel yet', flush=True)
        return
    young = [r for r in rows if r[1].startswith('young')]
    adult = [r for r in rows if not r[1].startswith('young')]
    order = young + adult

    fig, axes = plt.subplots(1, 3, figsize=(17, 4.8))

    # what each denominator does with age, in the cortex, per brain
    ax = axes[0]
    for i, r in enumerate(order):
        col = '#c0392b' if r[1].startswith('young') else '#555555'
        ax.plot(i, r[3], 'o', color=col, mfc='none', label='auto' if i == 0 else None)
        ax.plot(i, r[4], 's', color=col, label='SEP' if i == 0 else None)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([r[0].split('_')[0] for r in order], rotation=70, fontsize=7)
    ax.axvline(len(young) - 0.5, color='k', lw=0.6, ls=':')
    ax.set_ylabel('isocortex mean, raw counts'); ax.legend(fontsize=8)
    ax.set_title('the two denominators (open = auto, filled = SEP)', fontsize=10)
    ax.grid(lw=0.3, alpha=0.6)

    # how much of the nano difference each one would absorb
    ax = axes[1]
    for grp, col, lbl in ((young, '#c0392b', 'young'), (adult, '#555555', 'adult')):
        ax.plot([r[2] for r in grp], [r[4] for r in grp], 'o', color=col, label=lbl)
    ax.set_xlabel('isocortex mean, nano'); ax.set_ylabel('isocortex mean, SEP')
    ax.set_title('nano against SEP across brains', fontsize=10)
    ax.legend(fontsize=8); ax.grid(lw=0.3, alpha=0.6)

    # the part that matters: does the young-adult difference survive the swap
    ax = axes[2]
    stats = os.path.join(DATA, 'comparisons_v2', 'young_vs_adult', 'region_stats.csv')
    pairs = {}
    if os.path.exists(stats):
        for r in csv.DictReader(open(stats, encoding='utf-8')):
            if r['reading'] in ('ratio', 'sepratio'):
                pairs.setdefault(r['acronym'], {})[r['reading']] = float(r['diff_log2'])
    keys = [k for k, v in pairs.items() if len(v) == 2]
    if keys:
        x = np.array([pairs[k]['ratio'] for k in keys])
        y = np.array([pairs[k]['sepratio'] for k in keys])
        lim = float(max(np.abs(np.concatenate([x, y])).max(), 0.1)) * 1.1
        ax.plot([-lim, lim], [-lim, lim], '-', color='#bbbbbb', lw=1)
        ax.axhline(0, color='#dddddd', lw=0.8); ax.axvline(0, color='#dddddd', lw=0.8)
        ax.plot(x, y, 'o', ms=4, color='#c0392b')
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
        ax.set_title(f'young - adult, {len(keys)} structures (r = {np.corrcoef(x, y)[0, 1]:.2f})', fontsize=10)
    else:
        ax.text(0.5, 0.5, 'run v2_region_plot.py first', ha='center', va='center', transform=ax.transAxes)
        ax.set_title('young - adult per structure', fontsize=10)
    ax.set_xlabel('log2 difference, nano / auto'); ax.set_ylabel('log2 difference, nano / SEP')
    ax.grid(lw=0.3, alpha=0.6)

    fig.suptitle('Choosing the reference channel: autofluorescence measures tissue, SEP measures the receptor itself.\n'
                 'Points off the identity line on the right are structures whose answer depends on that choice', fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    save_figure(fig, os.path.join(OUT, '09_denominators.png'), dpi=95); plt.close(fig)


def write_index():
    with open(os.path.join(OUT, 'README.md'), 'w', encoding='utf-8') as fh:
        fh.write("""# processing diagnostics

One sheet per question, so every step of the v2 route can be checked by eye
rather than trusted. Regenerate with `v2_diagnostics.py`.

| sheet | what to look for | what would be wrong |
|---|---|---|
| `01_tissue_<mouse>.png` | red contour hugs the tissue edge; ventricles and the space around fibre tracts excluded; blue atlas outline roughly matches | mask eating into cortex, or spilling into the black surround |
| `02_levels_<mouse>.png` | off-tissue and tissue distributions separate at the marked threshold; cortex centred on 1.0 after scaling | the two distributions overlapping at the threshold, or a cortex peak far from 1 |
| `03_coverage.png` | each brain covers a contiguous run of planes; the group overlaps in the middle | a brain with holes, or a cohort where only one brain covers a whole region |
| `04_warp_<mouse>.png` | before and after look like the same brain; region means sit on the identity line | points off the line, i.e. the transform moving signal between regions |
| `05_cohort_n.png` | the n map is flat across most of the brain | large areas resting on one or two brains |
| `06_scaling.png` | backgrounds similar across brains; cortex means spread widely; the two channels track each other | a background far from the others (a brain imaged differently) |
| `07_route_agreement.png` | points on the identity line | a systematic offset, i.e. the warp biasing the comparison |
| `08_mask_vs_p6bis.png` | red and blue contours on top of each other | the new mask cutting into tissue, or reaching into the surround, where the old one does not |
| `09_denominators.png` | the two reference channels behave alike across brains; structures sit on the identity line | a structure whose young-adult difference flips sign with the denominator -- that result belongs to the denominator, not to the biology |
| `sep_channel/<mouse>.png` | written by `P4bis_add_sep_channel.m`: every slice recovered at r = 1.00000, and the registered DAPI of that run on top of the one P4 wrote | a slice below r = 0.99, or a residual shift above a fraction of a pixel: the SEP channel would not be sitting where NANO and AUTO sit |

The numbers behind these are in `../young_vs_adult/region_stats.csv` and
`../README.md`.
""")


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    ISO = set()
    with open(CSV_MAP, newline='', encoding='utf-8') as fh:
        for row in csv.DictReader(fh):
            if row['parcellation_term_set_name'] == 'division' and row['parcellation_term_acronym'] == 'Isocortex':
                ISO.add(int(row['parcellation_index']))
    ISO = list(ISO)
    ANN = {k: annotation_20(k) for k in {v[1] for v in MICE.values()}}
    mice = sys.argv[1:] or list(MICE)
    for mouse in mice:
        ann = ANN[MICE[mouse][1]]
        z = np.load(os.path.join(PER_MOUSE, mouse + '.npz'))
        sheet_tissue(mouse, ann, z)
        sheet_levels(mouse, ann, z)
        if MICE[mouse][1] != 'ccf':
            sheet_warp(mouse)
        print(f'{mouse:20s} sheets written', flush=True)
    if not sys.argv[1:]:
        sheet_coverage(); sheet_cohort_n(); sheet_scaling(); sheet_route_agreement()
        sheet_mask_vs_p6bis(); sheet_denominators(); write_index()
        print('cohort-level sheets and index written')
