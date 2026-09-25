"""
Does the nano map track receptor ON THE MEMBRANE, or just receptor?

This is the question the gene ranking alone cannot answer. Cacng8 outranking
Gria1 is suggestive, but Gria1, Cacng8, Dlg2 and Grip1 are all postsynaptic
genes, so no grouping of genes -- ours, GO's, or SynGO's -- distinguishes
"surface pool" from "total receptor". Both hypotheses predict a postsynaptic
map. The distinction lives in the CHANNELS, and v2_adult_arms.py puts the three
of them on one footing:

  sepauto    SEP / auto     total receptor, wherever it sits
  ratio      nano / auto    receptor at the surface
  sepratio   nano / SEP     the surface fraction, with auto divided out

Two tests are run on them, both stated before the numbers were computed.

**Test 1, the ordering across arms.** If SEP reports total receptor and nano
reports the surface pool, then against Gria1 expression the arms should order
sepauto > ratio > sepratio, and against trafficking and anchoring machinery the
order should reverse. Reported as a per-gene table and as the difference
sepratio - sepauto, which is the quantity the whole argument turns on.

**Test 2, partial correlation given Gria1.** Even a perfect surface map is
correlated with Gria1, because you cannot have surface receptor where there is
no receptor. So the sharper question is what is left after Gria1 is taken out:
rho(arm, gene | Gria1), on ranks. If the nano map were only reporting receptor
abundance, nothing should survive; if it reports the surface pool, the
machinery genes should still predict the residual, and should do so most in the
surface-fraction arm.

**THE PREMISE FAILED, AND THIS IS THE RECORD OF HOW.** Run first, then checked:
`SEP / auto` correlates with Gria1 at -0.11, which sent us to the channels
themselves. v2_sep_channel_check.py settles it -- in this fixed, cleared tissue
the green channel is mostly autofluorescence (rho 0.79 +- 0.04 against the
autofluorescence channel across all ten adults, against 0.26 for nano; dynamic
range 0.95 log2 against autofluo's 1.07 and nano's 1.93). So `sepauto` is not
total receptor, and `sepratio` is not a surface fraction -- it is nano divided
by a second autofluorescence-like channel, which is exactly why it tracks
`ratio` at rho 0.89 to 0.97 within every mouse.

Test 1 therefore measures nothing about membrane versus total, and is kept only
because it is how the failure was found. **Test 2 survives**, because it does
not use the green channel at all: Gria1 expression stands in for total receptor,
and the question is what the nano map still explains once it is removed. That is
a weaker reference than a real total-protein channel would be -- mRNA is not
protein -- but it is the one we have.

**What else can defeat this, and is not controlled here.**
 1. Gria1 mRNA as a stand-in for total receptor. Transport, translation and
    turnover all vary regionally, so a residual that trafficking genes predict
    is consistent with a surface map and also with a translation gradient.
    Only a real total-receptor channel separates those.
 2. Gria1's own ISH experiment is one mouse with its own section quality. A
    noisy Gria1 map weakens every correlation with it and inflates whatever
    survives the partial. Its coverage is printed for that reason.
 3. Spatial autocorrelation. Everything here is descriptive: the p-values that
    would go with these rho values are inflated (Fulcher 2021), which is why
    the arms are compared against EACH OTHER on the same genes and the same
    structures -- a paired contrast where that inflation largely cancels --
    rather than each being tested against zero.

Outputs, under data\\adult_v2\\arms:
  arm_gene_correlations.csv     arm x gene, plain and partial-on-Gria1
  arms_vs_genes.png             the two tests, drawn

  D:\\sep_histology\\code\\tools\\venv_atlas\\Scripts\\python.exe v2_ish_arms.py
"""

import csv
import os
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import rankdata, spearmanr, wilcoxon

DATA = r'D:\sep_histology\data'
ARMS_CSV = os.path.join(DATA, 'adult_v2', 'arms', 'region_means_arms.csv')
GENES = os.path.join(DATA, 'adult_v2', 'ish', 'gene_region_table.csv')
OUT = os.path.join(DATA, 'adult_v2', 'arms')

ARMS = ('sepauto', 'ratio', 'sepratio')
LABEL = {'sepauto': 'SEP / auto\ntotal receptor',
         'ratio': 'nano / auto\nsurface receptor',
         'sepratio': 'nano / SEP\nsurface fraction'}
CONTROL = 'Gria1'
MIN_ISH_VOXELS = 10
MIN_STRUCTURES = 50
# the genes the prediction is about: AMPAR anchoring and trafficking, named in
# gene_targets.csv, not chosen after seeing this result
MACHINERY = ('auxiliary', 'trafficking', 'scaffold')


def arm_profiles():
    """{arm: {structure: mean over the ten adults}}."""
    per = defaultdict(lambda: defaultdict(list))
    with open(ARMS_CSV, newline='', encoding='utf-8') as fh:
        for r in csv.DictReader(fh):
            per[r['arm']][r['structure']].append(float(r['log2_value']))
    return {arm: {s: float(np.mean(v)) for s, v in d.items()} for arm, d in per.items()}


def gene_profiles():
    """{gene: {structure: expression}}, {gene: category}, {gene: n structures}."""
    out, cat = defaultdict(dict), {}
    with open(GENES, newline='', encoding='utf-8') as fh:
        for r in csv.DictReader(fh):
            if int(r['n_voxels']) >= MIN_ISH_VOXELS:
                out[r['symbol']][r['structure']] = float(r['ish_mean'])
                cat[r['symbol']] = r['category']
    return out, cat


def partial_spearman(x, y, z):
    """Spearman of x and y with z removed, i.e. Pearson on the rank residuals."""
    rx, ry, rz = (rankdata(v).astype(float) for v in (x, y, z))
    zc = np.column_stack([rz, np.ones_like(rz)])
    resid = lambda r: r - zc @ np.linalg.lstsq(zc, r, rcond=None)[0]
    a, b = resid(rx), resid(ry)
    if a.std() == 0 or b.std() == 0:
        return np.nan
    return float(np.corrcoef(a, b)[0, 1])


def correlate(arms, genes, category, control):
    """One row per (arm, gene): rho, and rho with the control gene removed.

    Every arm is restricted to the SAME structures -- those measured in all
    three arms, in this gene, and in the control gene -- so the three numbers on
    a row are a paired comparison and not three different samples.
    """
    shared_arms = set.intersection(*[set(arms[a]) for a in ARMS])
    rows = []
    for gene, expr in sorted(genes.items()):
        if gene == control:
            common = sorted(shared_arms & set(expr))
        else:
            common = sorted(shared_arms & set(expr) & set(genes[control]))
        if len(common) < MIN_STRUCTURES:
            continue
        y = np.array([expr[s] for s in common])
        c = np.array([genes[control][s] for s in common])
        for arm in ARMS:
            x = np.array([arms[arm][s] for s in common])
            rho, _ = spearmanr(x, y)
            rows.append(dict(arm=arm, symbol=gene, category=category[gene],
                             n_structures=len(common), rho=float(rho),
                             rho_partial=(np.nan if gene == control
                                          else partial_spearman(x, y, c))))
    return rows


def by_gene(rows, field):
    """{gene: {arm: value}} for one column."""
    out = defaultdict(dict)
    for r in rows:
        out[r['symbol']][r['arm']] = r[field]
    return out


def report(rows, category):
    plain, partial = by_gene(rows, 'rho'), by_gene(rows, 'rho_partial')
    mach = sorted((g for g in plain if category[g] in MACHINERY),
                  key=lambda g: -plain[g]['sepratio'])

    print('\nTEST 1 -- the ordering across arms')
    print(f'  {"":12s} {"sepauto":>9s} {"ratio":>9s} {"sepratio":>9s}   '
          f'{"sepratio-sepauto":>17s}')
    named = [CONTROL] + mach[:6]
    for gene in named:
        p = plain[gene]
        print(f'  {gene:12s} {p["sepauto"]:+9.3f} {p["ratio"]:+9.3f} {p["sepratio"]:+9.3f}   '
              f'{p["sepratio"] - p["sepauto"]:+17.3f}')

    swing = {g: plain[g]['sepratio'] - plain[g]['sepauto'] for g in plain}
    m = [swing[g] for g in plain if category[g] in MACHINERY]
    print(f'\n  {CONTROL} swing towards the surface fraction: {swing[CONTROL]:+.3f}')
    print(f'  machinery genes (n = {len(m)}):  median {np.median(m):+.3f}, '
          f'{sum(1 for v in m if v > swing[CONTROL])} of {len(m)} above {CONTROL}')
    print('  the prediction is that Gria1 swings DOWN and the machinery swings up '
          'or falls less')

    print(f'\nTEST 2 -- with {CONTROL} partialled out, what is left')
    print(f'  {"":12s} {"sepauto":>9s} {"ratio":>9s} {"sepratio":>9s}')
    for gene in mach[:8]:
        p = partial[gene]
        print(f'  {gene:12s} {p["sepauto"]:+9.3f} {p["ratio"]:+9.3f} {p["sepratio"]:+9.3f}')
    for arm in ARMS:
        vals = [partial[g][arm] for g in plain if category[g] in MACHINERY
                and np.isfinite(partial[g][arm])]
        print(f'  {arm:12s} machinery median {np.median(vals):+.3f}, '
              f'{sum(1 for v in vals if v > 0)} of {len(vals)} positive')

    # the paired contrast: same genes, same structures, two arms
    a = np.array([plain[g]['sepratio'] for g in mach])
    b = np.array([plain[g]['sepauto'] for g in mach])
    stat, p = wilcoxon(a, b)
    print(f'\n  machinery, sepratio vs sepauto, paired over {len(mach)} genes: '
          f'median difference {np.median(a - b):+.3f}, Wilcoxon p = {p:.4f}')
    return plain, partial, mach


def figure(plain, partial, mach, category):
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.3))

    # left: every gene's rho across the three arms, Gria1 and machinery picked out
    ax = axes[0]
    x = np.arange(len(ARMS))
    for gene in sorted(plain):
        if gene == CONTROL or category[gene] in MACHINERY:
            continue
        ax.plot(x, [plain[gene][a] for a in ARMS], color='0.85', lw=0.7, zorder=1)
    for gene in mach:
        ax.plot(x, [plain[gene][a] for a in ARMS], color='#c0392b', lw=0.9,
                alpha=0.55, zorder=2)
    ax.plot(x, [plain[CONTROL][a] for a in ARMS], color='#1f3b73', lw=2.4, zorder=4,
            marker='o', ms=5)
    ax.annotate(CONTROL, (2, plain[CONTROL]['sepratio']), color='#1f3b73',
                fontsize=9, xytext=(6, -2), textcoords='offset points')
    ax.set_xticks(x); ax.set_xticklabels([LABEL[a] for a in ARMS], fontsize=7.5)
    ax.set_ylabel('Spearman with the arm, over structures', fontsize=8)
    ax.set_title('Test 1: every gene across the three arms\n'
                 'red = auxiliary / trafficking / scaffold', fontsize=9)

    # middle: the swing, machinery against everything else
    ax = axes[1]
    swing = {g: plain[g]['sepratio'] - plain[g]['sepauto'] for g in plain}
    groups = [[swing[g] for g in plain if category[g] not in MACHINERY and g != CONTROL],
              [swing[g] for g in mach]]
    rng = np.random.default_rng(0)
    for i, (vals, colour) in enumerate(zip(groups, ('0.65', '#c0392b'))):
        ax.scatter(np.full(len(vals), i) + rng.uniform(-.12, .12, len(vals)), vals,
                   s=15, facecolor=colour, edgecolor='0.25', linewidth=0.4, zorder=2)
        ax.plot([i - .28, i + .28], [np.median(vals)] * 2, color='0.15', lw=1.7, zorder=3)
    ax.axhline(swing[CONTROL], color='#1f3b73', lw=1.4, ls='--', zorder=1)
    ax.annotate(CONTROL, (1.35, swing[CONTROL]), color='#1f3b73', fontsize=8,
                va='bottom', ha='right')
    ax.axhline(0, color='0.8', lw=0.7, zorder=0)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f'other genes\n({len(groups[0])})', f'machinery\n({len(groups[1])})'],
                       fontsize=8)
    ax.set_ylabel('rho(surface fraction) - rho(total receptor)', fontsize=8)
    ax.set_title('Test 1: the swing towards the surface fraction', fontsize=9)

    # right: what survives partialling out Gria1
    ax = axes[2]
    for i, arm in enumerate(ARMS):
        vals = [partial[g][arm] for g in mach if np.isfinite(partial[g][arm])]
        ax.scatter(np.full(len(vals), i) + rng.uniform(-.12, .12, len(vals)), vals,
                   s=15, facecolor='#c0392b', edgecolor='0.25', linewidth=0.4, zorder=2)
        ax.plot([i - .28, i + .28], [np.median(vals)] * 2, color='0.15', lw=1.7, zorder=3)
    ax.axhline(0, color='0.8', lw=0.7, zorder=0)
    ax.set_xticks(range(len(ARMS)))
    ax.set_xticklabels([LABEL[a] for a in ARMS], fontsize=7.5)
    ax.set_ylabel(f'partial rho with {CONTROL} removed', fontsize=8)
    ax.set_title(f'Test 2: machinery genes, {CONTROL} partialled out', fontsize=9)

    for ax in axes:
        ax.tick_params(labelsize=7)
        for side in ('top', 'right'):
            ax.spines[side].set_visible(False)
    fig.suptitle('Membrane pool or total receptor?  The contrast is between channels, '
                 'not between genes.  Descriptive: arms are compared with each other on the '
                 'same genes and structures, never against zero.', fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    path = os.path.join(OUT, 'arms_vs_genes.png')
    fig.savefig(path, dpi=200); plt.close(fig)
    print(f'\n{path}')


def main():
    arms = arm_profiles()
    missing = [a for a in ARMS if a not in arms]
    if missing:
        raise SystemExit(f'arms missing from {ARMS_CSV}: {missing}. Run v2_adult_arms.py')
    genes, category = gene_profiles()
    if CONTROL not in genes:
        raise SystemExit(f'{CONTROL} is not in the gene table; it is the control here')
    print(f'{len(genes)} genes; {CONTROL} measured in {len(genes[CONTROL])} structures; '
          f'arms {ARMS}')

    rows = correlate(arms, genes, category, CONTROL)
    path = os.path.join(OUT, 'arm_gene_correlations.csv')
    with open(path, 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows({k: (f'{v:.4f}' if isinstance(v, float) else v) for k, v in r.items()}
                    for r in rows)
    print(f'{len(rows)} rows -> {path}')

    plain, partial, mach = report(rows, category)
    figure(plain, partial, mach, category)


if __name__ == '__main__':
    main()
