"""
What kind of gene sits at the top of the ranking? Annotation words, tested.

`v2_ish_compare.py` leaves one number per gene, and the only grouping available
so far is the `category` column we wrote by hand in gene_targets.csv while
choosing those genes. That column is coarse -- "trafficking" holds Cacng8 next
to Bsn and Syn1, which are presynaptic -- and using it to explain the ranking is
close to circular, because the categories are our own reasons for picking the
genes in the first place.

So this asks the same question with annotation nobody here wrote: every gene's
GO terms, from mygene.info, cached one JSON per gene. Two feature sets come out
of the same download:

  terms   whole GO terms, e.g. "regulation of receptor recycling"
          -- the grouping the hand-made categories were trying to be
  words   single words from those terms and from the gene's name, e.g.
          "recycling", "endosome", "postsynaptic", "channel"
          -- coarser, but it pools terms that say the same thing differently

For each feature held by at least MIN_GENES genes, the genes that have it are
compared against the genes that do not, on their correlation with the adult nano
map: Mann-Whitney on the rho values, the gap between the two medians as the
effect size with a percentile bootstrap interval around it, Benjamini-Hochberg
across features. Ranks rather than means because 95 genes is small and rho is
bounded. The interval matters more than it looks: sorted by effect size, a
five-gene group beats everything for free, and the whisker is what shows it.

**This is the null Fulcher 2021 measured at 875-fold inflation, and it is used
here with that understood.** Genes sharing a GO term are co-expressed, so their
rho values are not independent draws and the p-values are anticonservative by a
lot. The honest reading of the output is the EFFECT SIZE and the ordering, with
q as a tie-breaker between features of similar size -- not a claim that any one
term is significant. The figure says so on its face.

The second limit is the pool. These are 95 genes chosen for their relevance to
AMPA receptors, so "synaptic words win" is partly a description of how the panel
was built. Within-panel contrasts (which synaptic word beats which) are much
safer than any statement about synaptic genes in general; that one needs the
background panel in docs\\adult_ish_design.md.

Outputs, under data\\adult_v2\\ish:
  annotation\\<symbol>.json     the cached mygene.info record, one per gene
  feature_enrichment.csv       every feature, every reading
  ish_word_enrichment.png      the picture, for zref

  D:\\sep_histology\\code\\tools\\venv_atlas\\Scripts\\python.exe mapping\\v2_ish_words.py
"""

import csv
import json
import os
import re
import time
import urllib.parse
import urllib.request
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, false_discovery_control

from v2_paths import DATA
RHO = os.path.join(DATA, 'adult_v2', 'ish', 'gene_correlations.csv')
OUT = os.path.join(DATA, 'adult_v2', 'ish')
CACHE = os.path.join(OUT, 'annotation')

MYGENE = 'https://mygene.info/v3/query'
MIN_GENES = 5              # a feature needs this many genes to be tested
MAX_SHARE = 0.80           # and must not be on nearly every gene, which says nothing
MIN_WORD = 4               # characters; drops "of", "to", "ion" and friends
PLOT_READING = 'zref'
N_SHOWN = 12               # features per panel in the figure
N_BOOT = 2000              # resamples behind each feature's interval

# Words whose behaviour was predicted BEFORE this script existed, in the
# measurement-validation argument: if the nanobody reports surface GluA1 at
# glutamatergic postsynapses, the postsynaptic vocabulary should track the map
# and the presynaptic vocabulary should not -- and the panel carries both, 68
# genes against 49, so the split is not something the panel design forces.
PREDICTED = ('postsynaptic', 'presynaptic', 'axon', 'glutamatergic', 'dendritic',
             'spine', 'vesicle', 'inhibitory', 'gabaergic')
CONTRAST = ('postsynaptic', 'presynaptic')

# words that appear in GO terms as grammar rather than as meaning
STOP = {'process', 'activity', 'involved', 'response', 'positive', 'negative',
        'protein', 'cellular', 'from', 'into', 'with', 'that', 'this', 'other',
        'type', 'like', 'containing', 'complex'}


def fetch(symbol):
    """The mygene.info record for one mouse gene, cached on disk.

    Cached because it is the only step here that depends on a server being up,
    and because a re-run should not re-ask 95 times for an answer that does not
    change. Delete the folder to refresh.
    """
    os.makedirs(CACHE, exist_ok=True)
    path = os.path.join(CACHE, symbol + '.json')
    if os.path.exists(path):
        with open(path, encoding='utf-8') as fh:
            return json.load(fh)

    query = urllib.parse.urlencode({'q': f'symbol:{symbol}', 'species': 'mouse',
                                    'fields': 'symbol,name,summary,go', 'size': 5})
    with urllib.request.urlopen(f'{MYGENE}?{query}', timeout=60) as fh:
        hits = json.load(fh).get('hits', [])
    hits = [h for h in hits if h.get('symbol', '').lower() == symbol.lower()]
    record = hits[0] if hits else {}
    with open(path, 'w', encoding='utf-8') as fh:
        json.dump(record, fh, indent=1)
    time.sleep(0.2)                      # a courtesy to a free public API
    return record


def go_terms(record):
    """Every GO term on the record, as a set of strings.

    mygene returns a dict when a gene has one term in a branch and a list when
    it has several, so both shapes are handled. Evidence codes are kept in the
    cached JSON but not filtered on: dropping IEA would remove most of the
    annotation for the less-studied genes and bias the comparison towards the
    ones somebody has studied by hand.
    """
    terms = set()
    for branch in record.get('go', {}).values():
        for item in (branch if isinstance(branch, list) else [branch]):
            if isinstance(item, dict) and item.get('term'):
                terms.add(item['term'].strip().lower())
    return terms


def words_of(terms, name):
    """Single words from the GO terms and the gene's full name."""
    out = set()
    for text in list(terms) + [name or '']:
        for w in re.split(r'[^a-z]+', text.lower()):
            if len(w) < MIN_WORD or w in STOP:
                continue
            # a very light plural rule: enough to pool "spine"/"spines", and
            # held back from "across", "synapsis", "exocytosis" and the like
            if w.endswith('s') and len(w) > MIN_WORD + 1 and w[-2:] not in ('ss', 'us', 'is', 'as'):
                w = w[:-1]
            out.add(w)
    return out


def load_rho():
    """{reading: {gene: rho}} from what v2_ish_compare wrote."""
    per = defaultdict(dict)
    with open(RHO, newline='', encoding='utf-8') as fh:
        for r in csv.DictReader(fh):
            per[r['reading']][r['symbol']] = float(r['rho'])
    return per


def gap_interval(a, b, rng):
    """A percentile bootstrap interval for the difference of the two medians.

    Without it the panel is misleading: a five-gene group reaches a large median
    gap for nothing, and sorted by effect size those groups sit at the top. The
    interval puts that on the figure instead of in a footnote.
    """
    da = np.median(rng.choice(a, (N_BOOT, len(a))), axis=1)
    db = np.median(rng.choice(b, (N_BOOT, len(b))), axis=1)
    lo, hi = np.percentile(da - db, [2.5, 97.5])
    return float(lo), float(hi)


def test_features(features, rho):
    """One row per feature: the two medians, the gap, Mann-Whitney p, BH q.

    `features` is {feature: set of genes}. Genes without a rho for this reading
    are simply absent from both sides.
    """
    genes = set(rho)
    rng = np.random.default_rng(0)        # fixed, so a re-run gives the same intervals
    rows = []
    for feature, carriers in features.items():
        inside = sorted(carriers & genes)
        outside = sorted(genes - carriers)
        if not MIN_GENES <= len(inside) <= MAX_SHARE * len(genes) or len(outside) < MIN_GENES:
            continue
        a = np.array([rho[g] for g in inside])
        b = np.array([rho[g] for g in outside])
        _, p = mannwhitneyu(a, b, alternative='two-sided')
        lo, hi = gap_interval(a, b, rng)
        rows.append(dict(feature=feature, n_genes=len(inside),
                         median_with=float(np.median(a)),
                         median_without=float(np.median(b)),
                         gap=float(np.median(a) - np.median(b)),
                         gap_lo=lo, gap_hi=hi, p=float(p),
                         genes=' '.join(inside)))
    if rows:
        q = false_discovery_control([r['p'] for r in rows])
        for r, qi in zip(rows, q):
            r['q'] = float(qi)
    return sorted(rows, key=lambda r: -r['gap'])


def bars(ax, rows, title, xlim):
    """Top features by effect size; bar darkness is the q value, as elsewhere."""
    sel = rows[:N_SHOWN][::-1]
    y = np.arange(len(sel))
    shade = [str(np.clip(0.75 * min(r['q'], 1.0) ** 0.5, 0.05, 0.8)) for r in sel]
    gaps = np.array([r['gap'] for r in sel])
    err = np.abs(np.array([[r['gap'] - r['gap_lo'] for r in sel],
                           [r['gap_hi'] - r['gap'] for r in sel]]))
    ax.barh(y, gaps, color=shade, edgecolor='0.25', linewidth=0.5)
    ax.errorbar(gaps, y, xerr=err, fmt='none', ecolor='0.2', elinewidth=0.8, capsize=2)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r['feature'][:46]}  ({r['n_genes']})" for r in sel], fontsize=7)
    ax.axvline(0, color='0.3', lw=0.6)
    ax.set_xlim(xlim)
    ax.set_xlabel('median rho with the feature  -  median rho without', fontsize=8)
    ax.set_title(title, fontsize=9)
    ax.tick_params(axis='x', labelsize=7)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)


def strip(ax, rows, rho):
    """The predicted contrast, drawn as the distributions behind it.

    Deliberately not "whichever feature came top": that one is usually a
    five-gene group, and picking it after the fact is how a ranked list turns
    into a claim. These two words were named in advance.
    """
    by = {r['feature']: r for r in rows}
    rng = np.random.default_rng(0)
    labels, i = [], 0
    for word in CONTRAST:
        r = by.get(word)
        if r is None:
            continue
        carriers = set(r['genes'].split())
        for has, colour in ((False, '0.7'), (True, '#c0392b')):
            values = [rho[g] for g in sorted(carriers if has else set(rho) - carriers)]
            ax.scatter(np.full(len(values), i) + rng.uniform(-0.11, 0.11, len(values)),
                       values, s=13, facecolor=colour, edgecolor='0.25',
                       linewidth=0.4, zorder=2)
            ax.plot([i - 0.28, i + 0.28], [np.median(values)] * 2,
                    color='0.15', lw=1.6, zorder=3)
            labels.append(f'{"yes" if has else "no"}\n{len(values)}')
            i += 1
        ax.text(i - 1.5, ax.get_ylim()[0], f'{word}\ngap {r["gap"]:+.3f}, q = {r["q"]:.3f}',
                ha='center', va='bottom', fontsize=7.5)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=7.5)
    ax.set_ylabel(f'rho with the adult nano map ({PLOT_READING})', fontsize=8)
    ax.set_title('the contrast predicted in advance', fontsize=9)
    ax.tick_params(axis='y', labelsize=7)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)


def figure(terms, words, rho):
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.6),
                             gridspec_kw=dict(width_ratios=[1.3, 1.0, 0.8]))
    span = max([r['gap'] for r in terms[:N_SHOWN] + words[:N_SHOWN]] + [0.05])
    xlim = (min(-0.02, min(r['gap_lo'] for r in terms[:N_SHOWN] + words[:N_SHOWN]) - 0.02),
            max(r['gap_hi'] for r in terms[:N_SHOWN] + words[:N_SHOWN]) + 0.02)
    bars(axes[0], terms, 'GO terms', xlim)
    bars(axes[1], words, 'words, from GO terms and gene names', xlim)
    strip(axes[2], words, rho)
    fig.suptitle('What the top of the gene ranking is made of.  Bars are ordered by effect '
                 'size, whiskers are a 95% bootstrap interval, darkness is the BH q.'
                 '\nThe p behind that q is anticonservative, because genes sharing a term are '
                 'co-expressed (Fulcher 2021): read the size and the width, not the significance.',
                 fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    path = os.path.join(OUT, 'ish_word_enrichment.png')
    fig.savefig(path, dpi=200); plt.close(fig)
    print(f'\n{path}')


def main():
    rho = load_rho()
    genes = sorted(rho[PLOT_READING])
    print(f'{len(genes)} genes, readings {sorted(rho)}')

    term_genes, word_genes, no_go = defaultdict(set), defaultdict(set), []
    for gene in genes:
        record = fetch(gene)
        terms = go_terms(record)
        if not terms:
            no_go.append(gene)
        for t in terms:
            term_genes[t].add(gene)
        for w in words_of(terms, record.get('name', '')):
            word_genes[w].add(gene)
    print(f'{len(term_genes)} distinct GO terms, {len(word_genes)} distinct words'
          + (f'; no GO annotation for {no_go}' if no_go else ''))

    rows = []
    for reading in sorted(rho):
        for kind, features in (('term', term_genes), ('word', word_genes)):
            for r in test_features(features, rho[reading]):
                rows.append(dict(reading=reading, kind=kind, **r))

    path = os.path.join(OUT, 'feature_enrichment.csv')
    fields = ['reading', 'kind', 'feature', 'n_genes', 'median_with', 'median_without',
              'gap', 'gap_lo', 'gap_hi', 'p', 'q', 'genes']
    with open(path, 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows({k: (f'{v:.4g}' if isinstance(v, float) else v) for k, v in r.items()}
                    for r in rows)
    print(f'{len(rows)} rows -> {path}')

    terms = [r for r in rows if r['reading'] == PLOT_READING and r['kind'] == 'term']
    words = [r for r in rows if r['reading'] == PLOT_READING and r['kind'] == 'word']

    by = {r['feature']: r for r in words}
    print(f'\nthe words named in advance, for {PLOT_READING}')
    for word in PREDICTED:
        r = by.get(word)
        print(f'  {word:15s} ' + (f'gap {r["gap"]:+.3f}  n={r["n_genes"]:3d}  q={r["q"]:.3f}'
                                  if r else 'not tested (too few genes, or absent)'))
    for label, sel in (('GO terms', terms), ('words', words)):
        print(f'\ntop {label} for {PLOT_READING}  '
              f'(gap in median rho, n genes, BH q)')
        for r in sel[:N_SHOWN]:
            print(f'  {r["gap"]:+.3f} [{r["gap_lo"]:+.2f} {r["gap_hi"]:+.2f}]  '
                  f'n={r["n_genes"]:3d}  q={r["q"]:.3f}  {r["feature"]}')
    figure(terms, words, rho[PLOT_READING])


if __name__ == '__main__':
    main()
