"""
The adult nano map against every gene in the panel, structure by structure.

The question this exists to answer: **can zref replace the old pipeline's
normalisation?** The original route reached its adult map through a per-mouse
affine fit and a within-brain z-score; v2 reaches a comparable quantity with one
documented transform (see `zref` in v2_cohort). If the gene ranking comes out
the same, the elaborate version was not earning its keep and the method section
gets much shorter.

So this script does the simplest possible thing and reports how much it matters:

  nano side   one value per structure per adult, straight from
              young_vs_adult\\region_means_per_mouse.csv -- already computed,
              already audited -- averaged over the ten adults
  gene side   v2_ish_regions.py's table
  join        on structure name, which both sides key by
  statistic   Spearman over structures, per gene

It runs for every reading, not just zref, because they are free once the table
is loaded: `ratio` and `sepratio` (nano per unit autofluorescence and per unit
SEP), `cref`, `subref` and `zref`. Two of the three validation arms are here:

  nano      surface pool      -> ratio, cref, subref, zref
  nano/sep  surface fraction  -> sepratio: should look like trafficking
                                 machinery, and LESS like Gria1 than nano does

The third arm, SEP on its own (total receptor, predicted to look most like
Gria1), needs a reading that does not exist yet -- every current reading has
nano in the numerator. Adding `sep / auto` to v2_cohort would complete it; until
then the dissociation is measured as a shift between nano and nano/sep rather
than across all three.

**These are correlations, not tests.** Two brain maps agree partly because
everything is high in cortex and hippocampus; the conventional null for this is
inflated roughly 875-fold in mouse (Fulcher 2021, see docs\\adult_ish_design.md).
The ranking is descriptive and the output says so. A spatial null is the next
piece of work, not this one.

Outputs, all under data\\adult_v2\\ish:
  gene_correlations.csv          one row per gene and reading
  ish_old_vs_new.png             the diagnostic for the question above: every
                                 gene's old rho against its new one

  D:\\sep_histology\\code\\tools\\venv_atlas\\Scripts\\python.exe mapping\\v2_ish_compare.py
"""

import csv
import os
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

from v2_paths import DATA
NANO = os.path.join(DATA, 'comparisons_v2', 'young_vs_adult', 'region_means_per_mouse.csv')
GENES = os.path.join(DATA, 'adult_v2', 'ish', 'gene_region_table.csv')
OLD = os.path.join(DATA, 'comparisons', 'merged_naive_rws_vs_ish_summary_nosmooth',
                   'gene_panel_summary.csv')
OUT = os.path.join(DATA, 'adult_v2', 'ish')

ADULT_GROUPS = ('naive', 'rws')
READINGS = ('zref', 'cref', 'subref', 'ratio', 'sepratio')
MIN_ISH_VOXELS = 10        # a structure's gene value needs this many 200 um voxels
MIN_STRUCTURES = 50        # a gene needs this many shared structures to be correlated
MIN_GENES = 20             # and that many genes before the two rankings are compared

# the categories that carry the prediction, from gene_targets.csv
MACHINERY = ('auxiliary', 'trafficking', 'scaffold')


def adult_profile():
    """{reading: {structure: mean over the ten adults}} from the per-mouse table."""
    per = defaultdict(lambda: defaultdict(list))
    with open(NANO, newline='', encoding='utf-8') as fh:
        for r in csv.DictReader(fh):
            if r['group'] in ADULT_GROUPS:
                per[r['reading']][r['structure']].append(float(r['log2_value']))
    return {reading: {s: float(np.mean(v)) for s, v in d.items()}
            for reading, d in per.items()}


def gene_profiles():
    """{gene: {structure: expression}} and {gene: category}, for usable coverage."""
    out, cat = defaultdict(dict), {}
    with open(GENES, newline='', encoding='utf-8') as fh:
        for r in csv.DictReader(fh):
            if int(r['n_voxels']) >= MIN_ISH_VOXELS:
                out[r['symbol']][r['structure']] = float(r['ish_mean'])
                cat[r['symbol']] = r['category']
    return out, cat


def correlate(nano, genes, category):
    """One row per (reading, gene): Spearman over the structures they share."""
    rows = []
    for reading in READINGS:
        profile = nano[reading]
        for gene, expr in sorted(genes.items()):
            shared = sorted(set(profile) & set(expr))
            if len(shared) < MIN_STRUCTURES:
                continue
            rho, _ = spearmanr([profile[s] for s in shared], [expr[s] for s in shared])
            rows.append(dict(reading=reading, symbol=gene, category=category[gene],
                             rho=f'{rho:.4f}', n_structures=len(shared)))
    return rows


def old_ranking():
    """The MATLAB route's distance-weighted Spearman, {gene: rho}, if it is there."""
    if not os.path.exists(OLD):
        return {}
    with open(OLD, newline='', encoding='utf-8') as fh:
        return {r['symbol']: float(r['r_spearman_dw']) for r in csv.DictReader(fh)
                if r['r_spearman_dw'] not in ('', 'NaN')}


def rank_of(sel, gene):
    """(rank from the top, rho) for one gene in a sorted list of rows."""
    for i, r in enumerate(sel, 1):
        if r['symbol'] == gene:
            return i, float(r['rho'])
    return None, np.nan


def report(rows, old, category):
    """Everything the decision rests on, printed in the order it gets asked."""
    per_reading = {reading: sorted([r for r in rows if r['reading'] == reading],
                                   key=lambda r: -float(r['rho']))
                   for reading in READINGS}

    print('\ntop of each ranking, and where the two named genes sit')
    for reading, sel in per_reading.items():
        top = ', '.join(f"{r['symbol']} {float(r['rho']):+.2f}" for r in sel[:5])
        print(f'  {reading:9s} {top}')
        for gene in ('Cacng8', 'Gria1'):
            i, rho = rank_of(sel, gene)
            if i:
                print(f'  {"":9s}   {gene:8s} rank {i:3d}/{len(sel)}   rho {rho:+.3f}')

    if old:
        print('\ndoes the simple normalisation reproduce the old ranking?')
        for reading, sel in per_reading.items():
            new = {r['symbol']: float(r['rho']) for r in sel}
            both = sorted(set(old) & set(new))
            if len(both) < MIN_GENES:
                continue
            rho, _ = spearmanr([old[g] for g in both], [new[g] for g in both])
            print(f'  {reading:9s} rho = {rho:+.3f} between the two orderings, '
                  f'{len(both)} genes in common')

    # The claim made of the old route was not about a category average -- it was
    # that the genes topping the list are AMPAR anchoring and trafficking, with
    # Gria1 itself well down it. So that is what gets measured: how far Gria1 is
    # below the best machinery gene, and how many machinery genes are above it.
    print('\nthe dissociation, stated the way it was stated of the old route')
    print(f'  {"reading":9s} {"best machinery":>22s} {"Gria1":>15s} {"gap":>6s}  '
          f'machinery above Gria1')
    for reading, sel in per_reading.items():
        mach = [r for r in sel if r['category'] in MACHINERY]
        i, g = rank_of(sel, 'Gria1')
        if not mach or i is None:
            continue
        best, above = mach[0], sum(1 for r in mach if float(r['rho']) > g)
        print(f'  {reading:9s} {best["symbol"]:>12s} {float(best["rho"]):+8.3f} '
              f'{"rank " + str(i):>9s} {g:+6.3f} {float(best["rho"]) - g:6.3f}  '
              f'{above:>13d} / {len(mach)}')
    if old:
        sel = sorted(({'symbol': k, 'rho': f'{v:.4f}',
                       'category': category.get(k, '')} for k, v in old.items()),
                     key=lambda r: -float(r['rho']))
        i, g = rank_of(sel, 'Gria1')
        mach = [r for r in sel if r['category'] in MACHINERY]
        above = sum(1 for r in mach if float(r['rho']) > g)
        print(f'  {"old":9s} {mach[0]["symbol"]:>12s} {float(mach[0]["rho"]):+8.3f} '
              f'{"rank " + str(i):>9s} {g:+6.3f} {float(mach[0]["rho"]) - g:6.3f}  '
              f'{above:>13d} / {len(mach)}')


def figure(rows, old, category):
    """Old rho against new rho, gene by gene -- the picture of "did it change?"."""
    readings = [r for r in ('zref', 'cref', 'ratio', 'sepratio')
                if any(x['reading'] == r for x in rows)]
    fig, axes = plt.subplots(1, len(readings), figsize=(3.1 * len(readings), 3.4),
                             sharex=True, sharey=True)
    axes = np.atleast_1d(axes)

    everything = ([old[g] for g in old] +
                  [float(r['rho']) for r in rows if r['reading'] in readings])
    lim = [min(everything) - 0.06, max(everything) + 0.06]

    for ax, reading in zip(axes, readings):
        new = {r['symbol']: float(r['rho']) for r in rows if r['reading'] == reading}
        both = sorted(set(old) & set(new))
        x = np.array([old[g] for g in both])
        y = np.array([new[g] for g in both])
        mach = np.array([category.get(g, '') in MACHINERY for g in both])

        ax.axhline(0, color='0.85', lw=0.6, zorder=0)
        ax.axvline(0, color='0.85', lw=0.6, zorder=0)
        ax.plot(lim, lim, color='0.75', lw=0.8, ls='--', zorder=1)
        ax.scatter(x[~mach], y[~mach], s=14, facecolor='0.75', edgecolor='0.45',
                   linewidth=0.4, zorder=2)
        ax.scatter(x[mach], y[mach], s=20, facecolor='#c0392b', edgecolor='0.2',
                   linewidth=0.4, zorder=3)
        for gene, dx, dy in (('Cacng8', -34, 7), ('Gria1', 6, -3)):
            if gene in both:
                ax.annotate(gene, (old[gene], new[gene]), textcoords='offset points',
                            xytext=(dx, dy), fontsize=7)

        rho, _ = spearmanr(x, y)
        ax.set_title(f'{reading}\nrho = {rho:+.3f}, n = {len(both)}', fontsize=9)
        ax.set_xlabel('old route (affine + z-score)', fontsize=8)
        ax.set_xlim(lim); ax.set_ylim(lim)
        ax.tick_params(labelsize=7)
        for side in ('top', 'right'):
            ax.spines[side].set_visible(False)
    axes[0].set_ylabel('v2 route, one reading', fontsize=8)

    fig.suptitle('Gene ranking, old normalisation against v2   '
                 '(red = auxiliary, trafficking or scaffold; descriptive, not tested)',
                 fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    path = os.path.join(OUT, 'ish_old_vs_new.png')
    fig.savefig(path, dpi=200); plt.close(fig)
    print(f'\n{path}')


def main():
    nano = adult_profile()
    genes, category = gene_profiles()
    missing = [r for r in READINGS if r not in nano]
    if missing:
        raise SystemExit(f'readings missing from {NANO}: {missing}')
    print(f'{len(genes)} genes, {len(nano["zref"])} adult structures, '
          f'readings {list(READINGS)}')

    rows = correlate(nano, genes, category)
    if not rows:
        raise SystemExit('no gene shared enough structures with the nano table')

    path = os.path.join(OUT, 'gene_correlations.csv')
    with open(path, 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(f'{len(rows)} rows -> {path}')

    old = old_ranking()
    report(rows, old, category)
    if old:
        figure(rows, old, category)


if __name__ == '__main__':
    main()
