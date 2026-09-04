"""
Re-draw the P20-vs-adult figure and table from the cached volumes, without
redoing the five-minute transform. Same conventions as
compare_young_vs_adult_lrsum.py; see there for what the numbers mean.

Fixes over the first draft: panels drawn dorsal-up with the midline on the
right (they were transposed), and structures aggregated by name so a region
with several layers or parts in the ontology is one row, not several.

  D:\\sep_histology\\code\\tools\\venv_atlas\\Scripts\\python.exe replot_young_vs_adult_lrsum.py
"""

import csv
import os
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DATA = r'D:\sep_histology\data'
OUT = os.path.join(DATA, 'comparisons', 'young_P20_vs_adult_nano')
CSV_MAP = os.path.join(DATA, 'atlas', 'parcellation_to_parcellation_term_membership.csv')
CCF_AP0 = 180


def main():
    z = np.load(os.path.join(OUT, 'volumes_ccf20.npz'))
    adu_rel, p20_rel, log2r, both, ann20 = (z[k] for k in ('adu_rel', 'p20_rel', 'log2r', 'both', 'annot20'))
    inside = ann20 > 0
    eps = 1e-3

    # ------------------------------------------------ per-structure, by name
    names, acro, divi = {}, {}, {}
    with open(CSV_MAP, newline='', encoding='utf-8') as fh:
        for row in csv.DictReader(fh):
            idx = int(row['parcellation_index'])
            if row['parcellation_term_set_name'] == 'structure':
                names[idx] = row['parcellation_term_name']; acro[idx] = row['parcellation_term_acronym']
            elif row['parcellation_term_set_name'] == 'division':
                divi[idx] = row['parcellation_term_acronym']
    labs = ann20[both]
    sa, sp, cnt, meta = defaultdict(float), defaultdict(float), defaultdict(int), {}
    for lab, a, p in zip(labs, adu_rel[both], p20_rel[both]):
        key = names.get(int(lab), f'id{lab}')
        sa[key] += a; sp[key] += p; cnt[key] += 1
        meta[key] = (acro.get(int(lab), ''), divi.get(int(lab), ''))
    rows = []
    for key in cnt:
        if cnt[key] < 200:
            continue
        ma, mp = sa[key] / cnt[key], sp[key] / cnt[key]
        rows.append((key, meta[key][0], meta[key][1], ma, mp, np.log2((mp + eps) / (ma + eps)), cnt[key]))
    rows.sort(key=lambda r: -r[5])
    with open(os.path.join(OUT, 'region_log2ratio.csv'), 'w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(['structure', 'acronym', 'division', 'adult_rel', 'P20_rel', 'log2_P20_over_adult', 'voxels_20um'])
        for r in rows:
            w.writerow([r[0], r[1], r[2], f'{r[3]:.4f}', f'{r[4]:.4f}', f'{r[5]:.3f}', r[6]])
    print(f'{len(rows)} structures. Relative level = share of own isocortex median (adult n=10, P20 n=3).')
    cortex = sorted([r for r in rows if r[2] == 'Isocortex'], key=lambda r: -r[5])
    print('ISOCORTEX ranked by log2(P20/adult):')
    for r in cortex:
        tag = '  <-- visual' if r[1].startswith('VIS') else ('  <-- somatosensory' if r[1].startswith('SS') else '')
        print(f'  {r[1]:9s} {r[0][:34]:34s} {r[3]:6.2f} {r[4]:6.2f} {r[5]:+6.2f}{tag}')
    print(f'  {"acr":8s} {"structure":30s} {"adult":>6s} {"P20":>6s} {"log2":>6s}  division')
    print('MORE nanobody share in P20:')
    for r in rows[:18]:
        print(f'  {r[1]:8s} {r[0][:30]:30s} {r[3]:6.2f} {r[4]:6.2f} {r[5]:+6.2f}  {r[2]}')
    print('LESS nanobody share in P20:')
    for r in rows[-18:]:
        print(f'  {r[1]:8s} {r[0][:30]:30s} {r[3]:6.2f} {r[4]:6.2f} {r[5]:+6.2f}  {r[2]}')
    # the big structures, for orientation
    print('Largest structures:')
    for r in sorted(rows, key=lambda r: -r[6])[:12]:
        print(f'  {r[1]:8s} {r[0][:30]:30s} {r[3]:6.2f} {r[4]:6.2f} {r[5]:+6.2f}  {r[2]}')

    # ---------------------------------------------------------------- figure
    # planes: spread through the AP range where all 13 mice overlap (MG897 has
    # 25 sections, so that range is narrower than either cohort alone)
    cov = both.reshape(both.shape[0], -1).sum(1)
    ok = np.nonzero(cov > 0.5 * cov.max())[0]
    planes = [int(v) for v in np.linspace(ok[0], ok[-1], 5).round()]
    # colour scale on the cortex, not the hippocampus (which is >10x the cortex in both ages)
    vmax = 3.0
    fig, axes = plt.subplots(len(planes), 3, figsize=(13.5, 3.9 * len(planes)))
    for i, zc in enumerate(planes):
        panels = ((np.where(inside[zc], adu_rel[zc], np.nan), 'adult (naive + rws, n = 10)', 'hot', (0, vmax)),
                  (np.where(inside[zc], p20_rel[zc], np.nan), 'P20 carried into CCF (n = 3)', 'hot', (0, vmax)),
                  (log2r[zc], 'log2( P20 / adult ),  same share = 0', 'RdBu_r', (-2, 2)))
        for j, (im, title, cmap, lim) in enumerate(panels):
            ax = axes[i, j]
            # im is (DV, ML): dorsal at the top, lateral edge on the left, midline on the right
            h = ax.imshow(im, cmap=cmap, vmin=lim[0], vmax=lim[1], origin='upper',
                          interpolation='nearest', aspect='equal')
            ax.contour(inside[zc].astype(float), levels=[0.5], colors='w' if j < 2 else 'k', linewidths=0.35)
            ax.set_title(f'{title}   plane {zc * 2 + CCF_AP0} / 10 um', fontsize=9.5)
            ax.set_xticks([]); ax.set_yticks([])
            plt.colorbar(h, ax=ax, fraction=0.035, pad=0.01)
    fig.suptitle('Nano LR-sum, each cohort scaled to its own isocortex median; P20 carried into the adult CCF by the '
                 'DeMBA -> CCF transform.  Grey = outside both cohorts\' tissue.', fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'slices_adult_vs_P20.png'), dpi=105)
    print('figure rewritten')


if __name__ == '__main__':
    main()
