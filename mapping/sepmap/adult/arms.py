"""
The three arms of the measurement argument, as one region table per adult.

The gene ranking on its own cannot say whether the nanobody reports receptor ON
THE MEMBRANE or simply receptor. Gria1, Cacng8, Dlg2 and Grip1 are all
postsynaptic, so no grouping of genes -- ours, GO's or SynGO's -- separates the
two. What separates them is a contrast between CHANNELS, and the SEP channel
that run_add_sep_channel carried into registered space is the one that makes it possible:

  sepauto    SEP / autofluorescence     total receptor, however it is localised
  ratio      nano / autofluorescence    receptor at the surface
  sepratio   nano / SEP                 the surface FRACTION

The three share a denominator in a useful way -- in log space
`sepratio = ratio - sepauto` up to the difference between a mean of ratios and
a ratio of means -- so autofluorescence, and anything else common to both
channels, cancels out of the third. The predictions, stated here before the
numbers exist:

  sepauto  should track Gria1 expression best of the three
  ratio    should sit in between
  sepratio should track trafficking and anchoring machinery, and Gria1 least

**The premise, which the argument rested on, and which turned out to be wrong.**
SEP is superecliptic pHluorin on GluA1, so in a living cell it reports the
surface pool; this tissue is fixed, cleared and mounted, so its pH gradients are
gone and the green channel should report the receptor wherever it sits. It does
not. adult.sep_channel_check measures the green channel directly and finds it
tracks the autofluorescence channel at rho 0.79 +- 0.04 across all ten adults,
with a dynamic range of 0.95 log2 against nano's 1.93. Whatever tag survives the
protocol, autofluorescence dominates what is left.

The table below is therefore still correct arithmetic, and `sepauto` and
`sepratio` do not mean what their names promise: not total receptor and not a
surface fraction. They are kept because they are how that was established, and
because `ratio` and the consistency check are needed either way. Read
adult.sep_channel_check before using either of them for anything.

What this script produces is only the table: per adult, per structure, the three
arms in log2. The arithmetic is deliberately the same as young_vs_adult.region_plot's --
same 20 um annotation, same MIN_VOX, same mask-normalised smoothing of the
denominator, structures keyed by name over layer indices -- and it checks itself
against that script's output for the two arms they share. A silent divergence
there would invalidate any comparison downstream, so it is measured and printed
rather than assumed.

Outputs, under data\\adult_v2\\arms:
  region_means_arms.csv     arm x mouse x structure, log2
  arms_consistency.png      the self-check: against young_vs_adult.region_plot, and the
                            log-space identity between the three arms

  D:\\sep_histology\\code\\tools\\venv_atlas\\Scripts\\python.exe mapping\\run_adult_arms.py
"""

import csv
import math
import os
from collections import defaultdict

import numpy as np
from scipy.ndimage import gaussian_filter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from sepmap.volumes.per_mouse import annotation_20, MICE, CSV_MAP, OUT as PER_MOUSE, DATA
from sepmap.volumes.cohort import RATIO_CLIP, NAIVE, RWS

OUT = os.path.join(DATA, 'adult_v2', 'arms')
EXISTING = os.path.join(DATA, 'comparisons_v2', 'young_vs_adult', 'region_means_per_mouse.csv')

MIN_VOX = 250              # 20 um voxels, as in young_vs_adult.region_plot
ADULTS = NAIVE + RWS
ARMS = ('sepauto', 'ratio', 'sepratio')
LABEL = {'sepauto': 'SEP / autofluorescence   (total receptor)',
         'ratio': 'nano / autofluorescence   (surface receptor)',
         'sepratio': 'nano / SEP   (surface fraction)'}
# the two arms young_vs_adult.region_plot already computes, and must agree with
SHARED = {'ratio': 'ratio', 'sepratio': 'sepratio'}


def structure_meta():
    """parcellation_index -> (name, acronym, division acronym)."""
    names, acro, divi = {}, {}, {}
    with open(CSV_MAP, newline='', encoding='utf-8') as fh:
        for row in csv.DictReader(fh):
            idx = int(row['parcellation_index'])
            if row['parcellation_term_set_name'] == 'structure':
                names[idx] = row['parcellation_term_name']
                acro[idx] = row['parcellation_term_acronym']
            elif row['parcellation_term_set_name'] == 'division':
                divi[idx] = row['parcellation_term_acronym']
    return names, acro, divi


def per_unit(num, ref, tissue):
    """num / ref voxel by voxel, exactly as volumes.cohort and young_vs_adult.region_plot do it.

    The denominator is smoothed by one 20 um voxel so that a single dark voxel
    cannot blow the ratio up, and the smoothing is normalised by the mask so
    tissue at the edge is not divided by the black outside it.
    """
    ref_s = gaussian_filter(np.where(tissue, ref, 0), 1.0) / np.maximum(
        gaussian_filter(tissue.astype(np.float32), 1.0), 1e-3)
    r = np.where(ref_s > 0, num / np.maximum(ref_s, 1e-3), 0)
    return np.clip(r, -RATIO_CLIP, RATIO_CLIP)


def mouse_table(mouse, names):
    """{structure name: (n voxels, mean per arm)} for one adult."""
    z = np.load(os.path.join(PER_MOUSE, mouse + '.npz'))
    if 'sep' not in z.files:
        raise SystemExit(f'{mouse}: no SEP channel. Run run_add_sep_channel.m, '
                         f'then run_per_mouse.py, for this brain.')
    sig = z['sig'].astype(np.float32)
    auto = z['auto'].astype(np.float32)
    sep = z['sep'].astype(np.float32)
    tissue = z['tissue']
    ann = annotation_20(MICE[mouse][1])

    vol = {'sepauto': per_unit(sep, auto, tissue),
           'ratio': per_unit(sig, auto, tissue),
           'sepratio': per_unit(sig, sep, tissue)}

    lab = ann[tissue]
    nlab = int(ann.max()) + 1
    n = np.bincount(lab, minlength=nlab)
    sums = {arm: np.bincount(lab, weights=vol[arm][tissue], minlength=nlab) for arm in ARMS}

    acc = defaultdict(lambda: [0] + [0.0] * len(ARMS))
    for idx in np.nonzero(n)[0]:
        if idx == 0:
            continue
        a = acc[names.get(int(idx), f'id{idx}')]
        a[0] += int(n[idx])
        for j, arm in enumerate(ARMS, 1):
            a[j] += sums[arm][idx]
    return {k: (v[0], {arm: v[j] / v[0] for j, arm in enumerate(ARMS, 1)})
            for k, v in acc.items() if v[0] >= MIN_VOX}


def check_against_existing(rows):
    """The two arms young_vs_adult.region_plot also computes must come out identical.

    They are computed here from the same per-mouse files with the same
    arithmetic, so anything above rounding means the two scripts have drifted
    apart and nothing downstream can be trusted.
    """
    if not os.path.exists(EXISTING):
        print('no existing table to check against -- skipped')
        return {}
    ours = {(r['arm'], r['mouse'], r['structure']): float(r['log2_value']) for r in rows}
    diffs = defaultdict(list)
    with open(EXISTING, newline='', encoding='utf-8') as fh:
        for r in csv.DictReader(fh):
            key = (SHARED.get(r['reading'], ''), r['mouse'], r['structure'])
            if key in ours:
                diffs[key[0]].append(abs(ours[key] - float(r['log2_value'])))
    # that table is written to four decimals, so half of the last digit is the
    # most the two can differ by if they are computing the same thing
    tol = 0.5e-4
    print(f'\nagainst run_region_plot (bound is {tol:.1e}, half the last digit it stores):')
    ok = True
    for arm, d in sorted(diffs.items()):
        worst = max(d)
        ok &= worst <= tol
        print(f'  {arm:9s} n = {len(d):5d}   max |difference| = {worst:.3e}   '
              f'{"agrees" if worst <= tol else "DRIFTED"}')
    if not ok:
        raise SystemExit('the two scripts no longer compute the same thing -- stop here')
    return diffs


def figure(per, diffs):
    """The self-check, drawn: the agreement, and the log-space identity."""
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 3.9))

    ax = axes[0]
    if diffs:
        for i, (arm, d) in enumerate(sorted(diffs.items())):
            ax.scatter(np.full(len(d), i) + np.random.default_rng(0).uniform(-.12, .12, len(d)),
                       np.maximum(d, 1e-17), s=5, facecolor='0.4', edgecolor='none')
        ax.set_yscale('log'); ax.set_xticks(range(len(diffs)))
        ax.set_xticklabels(sorted(diffs), fontsize=8)
        ax.axhline(1e-9, color='#c0392b', lw=0.8, ls='--')
        ax.set_ylabel('|this script - run_region_plot|  (log2 units)', fontsize=8)
    ax.set_title('the two shared arms agree', fontsize=9)

    # log(nano/sep) should equal log(nano/auto) - log(sep/auto), except that
    # each arm is a mean of voxelwise ratios rather than a ratio of means
    gap = []
    for mouse, table in per.items():
        for k, (_, m) in table.items():
            if min(m.values()) > 0:
                gap.append(math.log2(m['sepratio'])
                           - (math.log2(m['ratio']) - math.log2(m['sepauto'])))
    axes[1].hist(gap, bins=60, color='0.6', edgecolor='0.3', linewidth=0.4)
    axes[1].axvline(0, color='#c0392b', lw=0.9)
    axes[1].set_xlabel('log2(nano/SEP)  -  [log2(nano/auto) - log2(SEP/auto)]', fontsize=8)
    axes[1].set_title(f'Jensen gap: median {np.median(gap):+.3f}, '
                      f'p5-p95 {np.percentile(gap, 5):+.2f} to {np.percentile(gap, 95):+.2f}',
                      fontsize=9)
    axes[1].set_ylabel('structures x mice', fontsize=8)

    # and the two channels against each other, one dot per structure
    ax = axes[2]
    mouse = sorted(per)[0]
    t = per[mouse]
    ks = [k for k in t if min(t[k][1].values()) > 0]
    ax.scatter([math.log2(t[k][1]['sepauto']) for k in ks],
               [math.log2(t[k][1]['ratio']) for k in ks],
               s=9, facecolor='0.55', edgecolor='0.25', linewidth=0.3)
    ax.set_xlabel('log2 SEP / auto', fontsize=8)
    ax.set_ylabel('log2 nano / auto', fontsize=8)
    ax.set_title(f'the two channels, {mouse}', fontsize=9)

    for ax in axes:
        ax.tick_params(labelsize=7)
        for side in ('top', 'right'):
            ax.spines[side].set_visible(False)
    fig.tight_layout()
    path = os.path.join(OUT, 'arms_consistency.png')
    fig.savefig(path, dpi=200); plt.close(fig)
    print(f'{path}')


def main():
    os.makedirs(OUT, exist_ok=True)
    names, acro, divi = structure_meta()
    meta = {nm: (acro[idx], divi.get(idx, '')) for idx, nm in names.items()}
    group = {m: ('naive' if m in NAIVE else 'rws') for m in ADULTS}

    per, rows = {}, []
    for mouse in ADULTS:
        per[mouse] = mouse_table(mouse, names)
        print(f'{mouse:20s} {len(per[mouse])} structures', flush=True)
        for k, (n, m) in per[mouse].items():
            for arm in ARMS:
                if m[arm] > 0:
                    rows.append(dict(arm=arm, group=group[mouse], mouse=mouse, structure=k,
                                     acronym=meta.get(k, ('', ''))[0],
                                     division=meta.get(k, ('', ''))[1],
                                     n_vox20=n, log2_value=f'{math.log2(m[arm]):.6f}'))

    path = os.path.join(OUT, 'region_means_arms.csv')
    with open(path, 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(f'\n{len(rows):,} rows -> {path}')

    diffs = check_against_existing(rows)
    figure(per, diffs)
