"""
How reliable is one Allen map? Measured, from the genes measured twice.

Every correlation in this analysis has been between the nano map and ONE in situ
experiment per gene -- one mouse, one set of sections, one staining run. A
correlation is capped by the reliability of both things being correlated, so a
gene whose map is noisy looks uninformative whether it is or not. That has been
an open caveat since the first version of this comparison ("check Gria1's own
experiment first") and it has never been a number.

The ontology panel makes it one. Of its genes, a large fraction have both a
coronal and a sagittal experiment: two independent passes at the same gene,
different animals, different sectioning. Correlating those two profiles is a
test-retest reliability for a single Allen map.

Two things come out of it:

  reliability   per gene, the Spearman between its experiments. This is what
                every correlation involving that gene is capped by, and it is
                also a quality filter that is measured instead of assumed.
  merged map    the average of a gene's experiments, which is a better estimate
                than either. Averaged on RANKS, because expression energy is
                not calibrated across experiments -- two runs of the same gene
                can differ by a constant factor for reasons that have nothing
                to do with the brain.

**Coronal against sagittal is the honest pairing and also the harshest.** The
two planes are reconstructed differently and a sagittal experiment covers one
hemisphere, so their agreement is a lower bound on what two experiments of the
same plane would give. Where a gene has two of the same plane, that is reported
separately, and the difference between the two is worth looking at before
leaning on either.

Outputs, under data\\adult_v2\\ish:
  gene_reliability.csv          per gene: n experiments, reliability, coverage
  gene_region_table_merged.csv  one profile per gene, replicates averaged
  ish_reliability.png           the distribution, and what it depends on

  D:\\sep_histology\\code\\tools\\venv_atlas\\Scripts\\python.exe mapping\\run_ish_reliability.py
"""

import csv
import itertools
import os
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import rankdata, spearmanr

from sepmap.config import DATA, SETTINGS
OUT = os.path.join(DATA, 'adult_v2', 'ish')
# The table to read is the one run_ish_regions wrote for one of the panel passes
# named in settings.toml ([ish_panels]); by default the ontology panel, whose
# genes include those with more than one experiment.
ISH_PANELS = SETTINGS['ish_panels']
DEFAULT_PANEL = 'ontology'

MIN_ISH_VOXELS = 10
MIN_SHARED = 50            # structures two experiments must share to be compared


def load(path):
    """{gene: {experiment: (plane, role, {structure: energy})}}."""
    per = defaultdict(lambda: defaultdict(dict))
    meta = {}
    with open(path, newline='', encoding='utf-8') as fh:
        for r in csv.DictReader(fh):
            if int(r['n_voxels']) < MIN_ISH_VOXELS:
                continue
            per[r['symbol']][r['experiment_id']][r['structure']] = float(r['ish_mean'])
            meta[(r['symbol'], r['experiment_id'])] = (r.get('plane', ''), r['category'])
    return per, meta


def pair_reliability(profiles):
    """Spearman between every pair of a gene's experiments, and the pairing."""
    out = []
    for a, b in itertools.combinations(sorted(profiles), 2):
        shared = sorted(set(profiles[a]) & set(profiles[b]))
        if len(shared) < MIN_SHARED:
            continue
        rho, _ = spearmanr([profiles[a][s] for s in shared],
                           [profiles[b][s] for s in shared])
        out.append((a, b, float(rho), len(shared)))
    return out


def merge(profiles):
    """One profile per gene: the mean of its experiments' ranks.

    Ranks rather than values because expression energy carries an arbitrary
    per-experiment scale. Structures are kept if at least one experiment has
    them, and each experiment's ranks are put on 0-1 first so that experiments
    covering different numbers of structures can still be averaged.
    """
    acc = defaultdict(list)
    for exp, prof in profiles.items():
        keys = sorted(prof)
        r = rankdata([prof[k] for k in keys])
        r = (r - 1) / max(len(r) - 1, 1)
        for k, v in zip(keys, r):
            acc[k].append(v)
    return {k: float(np.mean(v)) for k, v in acc.items()}, \
           {k: len(v) for k, v in acc.items()}


def main(panel_name=DEFAULT_PANEL):
    TABLE = ISH_PANELS[panel_name]['table']
    path = os.path.join(OUT, TABLE)
    if not os.path.exists(path):
        raise SystemExit(f'{path} not found -- run run_ish_regions.py with '
                         f'--panel {panel_name} first')
    per, meta = load(path)
    print(f'{len(per)} genes, {sum(len(v) for v in per.values())} experiments, '
          f'from {TABLE}')

    rows, merged_rows = [], []
    rel_by_pairing = defaultdict(list)
    for gene, profiles in sorted(per.items()):
        pairs = pair_reliability(profiles)
        rel = float(np.median([p[2] for p in pairs])) if pairs else np.nan
        planes = sorted({meta[(gene, e)][0] for e in profiles})
        role = meta[(gene, sorted(profiles)[0])][1]
        for a, b, rho, n in pairs:
            pa, pb = meta[(gene, a)][0], meta[(gene, b)][0]
            rel_by_pairing['same plane' if pa == pb else 'coronal vs sagittal'].append(rho)

        prof, n_exp = merge(profiles)
        level = float(np.median([v for p in profiles.values() for v in p.values()]))
        rows.append(dict(symbol=gene, role=role, n_experiments=len(profiles),
                         planes=' '.join(planes), n_pairs=len(pairs),
                         reliability=('' if np.isnan(rel) else f'{rel:.4f}'),
                         n_structures=len(prof), median_energy=f'{level:.4g}'))
        for structure, value in prof.items():
            merged_rows.append(dict(symbol=gene, role=role, structure=structure,
                                    rank_mean=f'{value:.6f}',
                                    n_experiments=n_exp[structure]))

    p1 = os.path.join(OUT, 'gene_reliability.csv')
    with open(p1, 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    p2 = os.path.join(OUT, 'gene_region_table_merged.csv')
    with open(p2, 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=list(merged_rows[0].keys()))
        w.writeheader(); w.writerows(merged_rows)
    print(f'{len(rows)} genes -> {p1}\n{len(merged_rows):,} rows -> {p2}')

    have = [r for r in rows if r['reliability'] != '']
    rel = np.array([float(r['reliability']) for r in have])
    print(f'\n{len(have)} genes measured more than once')
    print(f'  reliability  median {np.median(rel):+.3f}   '
          f'quartiles [{np.percentile(rel, 25):+.3f} {np.percentile(rel, 75):+.3f}]   '
          f'{(rel < 0.3).sum()} below 0.3')
    for pairing, v in sorted(rel_by_pairing.items()):
        print(f'  {pairing:20s} n={len(v):4d}  median {np.median(v):+.3f}')

    print('\n  the genes the test turns on:')
    for gene in ('Gria1', 'Gria2', 'Gria3', 'Gria4', 'Cacng8', 'Dlg4', 'Dlg2', 'Nsf'):
        r = next((x for x in rows if x['symbol'] == gene), None)
        if r:
            print(f'    {gene:8s} {r["n_experiments"]} experiments ({r["planes"]}), '
                  f'reliability {r["reliability"] or "n/a"}')

    figure(rows, rel, rel_by_pairing, per, meta)


def figure(rows, rel, rel_by_pairing, per, meta):
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.0))

    axes[0].hist(rel, bins=40, color='0.7', edgecolor='0.35', linewidth=0.4)
    axes[0].axvline(float(np.median(rel)), color='#c0392b', lw=1.6)
    axes[0].set_xlabel('Spearman between two experiments of the same gene', fontsize=8)
    axes[0].set_ylabel('genes', fontsize=8)
    axes[0].set_title(f'how reliable one Allen map is\nmedian {np.median(rel):+.2f}, '
                      f'n = {len(rel)}', fontsize=9)

    ax = axes[1]
    have = [r for r in rows if r['reliability'] != '']
    x = np.log10([float(r['median_energy']) + 1e-3 for r in have])
    y = [float(r['reliability']) for r in have]
    ax.scatter(x, y, s=10, facecolor='0.55', edgecolor='0.25', linewidth=0.3)
    ax.axhline(0, color='0.85', lw=0.7)
    ax.set_xlabel('log10 median expression energy', fontsize=8)
    ax.set_ylabel('reliability', fontsize=8)
    ax.set_title(f'a quiet gene is an unreliable one\nrho = '
                 f'{spearmanr(x, y).statistic:+.2f}', fontsize=9)

    ax = axes[2]
    rng = np.random.default_rng(0)
    keys = sorted(rel_by_pairing)
    for i, k in enumerate(keys):
        v = rel_by_pairing[k]
        ax.scatter(np.full(len(v), i) + rng.uniform(-.12, .12, len(v)), v, s=11,
                   facecolor='0.6', edgecolor='0.25', linewidth=0.3, zorder=2)
        ax.plot([i - .28, i + .28], [np.median(v)] * 2, color='#c0392b', lw=1.8, zorder=3)
    ax.set_xticks(range(len(keys)))
    ax.set_xticklabels([f'{k}\n({len(rel_by_pairing[k])})' for k in keys], fontsize=8)
    ax.axhline(0, color='0.85', lw=0.7)
    ax.set_ylabel('reliability', fontsize=8)
    ax.set_title('which pairing was available', fontsize=9)

    for ax in axes:
        ax.tick_params(labelsize=7)
        for side in ('top', 'right'):
            ax.spines[side].set_visible(False)
    fig.suptitle('Every correlation in this analysis is capped by these numbers: '
                 'a gene measured once, badly, cannot correlate with anything', fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    path = os.path.join(OUT, 'ish_reliability.png')
    fig.savefig(path, dpi=200); plt.close(fig)
    print(f'\n{path}')
