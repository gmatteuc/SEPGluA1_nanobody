"""
Is the nano map more than receptor abundance, and more than synaptic density?

These are the two things a sceptical reader will say the map is, and ranking
genes cannot answer either -- v2_ish_panel_test.py showed the ranking is flat
across postsynaptic gene classes, which is exactly what you would expect if the
map were only tracking how much synapse a region has. So the question is asked
the other way round here: **take both explanations out, and ask whether what is
left is signal or noise.**

    nano zref  ~  AMPA abundance  +  synaptic density        -> residual
    is the residual reproducible across independent animals?

That last line is the whole argument. A residual is always non-zero; the only
thing that makes it interesting is whether two halves of the cohort produce the
SAME residual. Noise does not survive that, and anything that does is real
spatial structure in surface GluA1 that neither explanation accounts for.

**The four covariates, and why each is there.**

  abundance   Gria1-4, the mean of their rank profiles. The lower bar: if the
              map is just "where AMPA receptor is expressed", this removes it.
  markers     the canonical synaptic markers in the panel -- Syp, Syn1, Vamp2,
              Bsn, Syt1 presynaptically, Dlg4, Homer1, Shank2, Shank3, Nlgn1,
              Camk2a postsynaptically. The higher bar, in its most legible form.
  psd_pc1     the first principal component of all 300 postsynaptic-density
              genes: the dominant axis of postsynaptic gene expression across
              the brain, taken from the data rather than named by us. A second,
              less arguable version of the same bar.
  autofluo    the autofluorescence channel of the SAME ten brains, per
              structure. This one is worth more than it looks: it is an
              internal measure of tissue and neuropil density, in the same
              sections, with no cross-modality registration between it and the
              nano channel. Every gene covariate is an Allen mouse warped into
              our atlas; this one is the actual tissue the nano signal came from.

**The ceiling matters as much as the model.** A map measured in ten animals has
its own reliability, and no covariate can explain variance that is not there.
So the cohort is split in half every possible way (126 splits of 10 into 5 and
5), the two half-maps correlated, and Spearman-Brown applied: that is the most
any predictor could reach. Explained variance is then reported as a fraction of
the RELIABLE variance, which is the honest denominator.

**What this can and cannot show.** It can show that the nano map contains
reproducible spatial structure that receptor-abundance and synaptic-density
proxies do not account for -- which is the interesting claim, and the one a
reader will want. It cannot show that the leftover IS the surface fraction:
mRNA is an imperfect proxy for protein, and a residual is only ever "not
explained by the things we put in". Anything that varies across the brain and
was not in the model lands in it.

Outputs, under data\\adult_v2\\beyond:
  variance_partition.csv     what each covariate set explains, against the ceiling
  residual_by_structure.csv  where the map exceeds and falls short of prediction
  beyond_density.png         the ceiling, the partition, and the reproducibility

  D:\\sep_histology\\code\\tools\\venv_atlas\\Scripts\\python.exe v2_beyond_density.py
"""

import csv
import itertools
import math
import os
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import rankdata, spearmanr

from v2_per_mouse import annotation_20, MICE, CSV_MAP, OUT as PER_MOUSE, DATA
from v2_cohort import NAIVE, RWS

NANO = os.path.join(DATA, 'comparisons_v2', 'young_vs_adult', 'region_means_per_mouse.csv')
ISH = os.path.join(DATA, 'adult_v2', 'ish', 'gene_region_table_merged.csv')
OUT = os.path.join(DATA, 'adult_v2', 'beyond')

ADULTS = NAIVE + RWS
READING = 'zref'
MIN_VOX = 250
HALF = 5

PER_EXPERIMENT = os.path.join(DATA, 'adult_v2', 'ish', 'gene_region_table_panel.csv')

SUBUNITS = ('Gria1', 'Gria2', 'Gria3', 'Gria4')
MARKERS = ('Syp', 'Syn1', 'Vamp2', 'Bsn', 'Syt1',            # presynaptic
           'Dlg4', 'Homer1', 'Shank2', 'Shank3', 'Nlgn1', 'Camk2a')   # postsynaptic


def nano_per_mouse():
    """{mouse: {structure: zref}} for the ten adults."""
    per = defaultdict(dict)
    with open(NANO, newline='', encoding='utf-8') as fh:
        for r in csv.DictReader(fh):
            if r['reading'] == READING and r['mouse'] in ADULTS:
                per[r['mouse']][r['structure']] = float(r['log2_value'])
    return per


def auto_per_mouse(structures):
    """{mouse: {structure: log2 mean autofluorescence}}, from the same brains."""
    names = {}
    with open(CSV_MAP, newline='', encoding='utf-8') as fh:
        for row in csv.DictReader(fh):
            if row['parcellation_term_set_name'] == 'structure':
                names[int(row['parcellation_index'])] = row['parcellation_term_name']
    per = {}
    for mouse in ADULTS:
        z = np.load(os.path.join(PER_MOUSE, mouse + '.npz'))
        tissue = z['tissue']
        ann = annotation_20(MICE[mouse][1])
        lab = ann[tissue]
        n = np.bincount(lab, minlength=int(ann.max()) + 1)
        s = np.bincount(lab, weights=z['auto'].astype(np.float32)[tissue], minlength=len(n))
        acc = defaultdict(lambda: [0, 0.0])
        for i in np.nonzero(n)[0]:
            if i == 0:
                continue
            a = acc[names.get(int(i), f'id{i}')]
            a[0] += int(n[i]); a[1] += s[i]
        per[mouse] = {k: math.log2(t / c) for k, (c, t) in acc.items()
                      if c >= MIN_VOX and t > 0 and k in structures}
        print(f'  {mouse:20s} {len(per[mouse])} structures', flush=True)
    return per


def gene_profiles():
    per, role = defaultdict(dict), {}
    with open(ISH, newline='', encoding='utf-8') as fh:
        for r in csv.DictReader(fh):
            per[r['symbol']][r['structure']] = float(r['rank_mean'])
            role[r['symbol']] = r['role']
    return per, role


def composite(genes, expr, structures):
    return np.mean([rankdata([expr[g][s] for s in structures]) for g in genes], axis=0)


def first_pc(genes, expr, structures):
    """The dominant axis of a gene set, on ranks, sign-fixed to the set mean."""
    m = np.array([rankdata([expr[g][s] for s in structures]) for g in genes])
    m = (m - m.mean(axis=1, keepdims=True)) / m.std(axis=1, keepdims=True)
    u, s, vt = np.linalg.svd(m - m.mean(axis=0), full_matrices=False)
    pc = vt[0]
    if np.corrcoef(pc, m.mean(axis=0))[0, 1] < 0:
        pc = -pc
    return pc, float(s[0] ** 2 / (s ** 2).sum())


def covariate_reliability(genes, structures):
    """Test-retest of a COMPOSITE, by building it twice from different experiments.

    The obvious objection to a large residual is that the covariates are noisy:
    a predictor measured badly cannot remove the variance it ought to. This
    answers it with the same trick as v2_ish_reliability -- genes measured more
    than once are split, one experiment into each of two composites, and the two
    composites correlated. The number it returns is the ceiling on what this
    covariate could explain even if it were measured perfectly.
    """
    runs = defaultdict(lambda: defaultdict(dict))
    with open(PER_EXPERIMENT, newline='', encoding='utf-8') as fh:
        for r in csv.DictReader(fh):
            if r['symbol'] in genes and int(r['n_voxels']) >= 10:
                runs[r['symbol']][r['experiment_id']][r['structure']] = float(r['ish_mean'])

    chosen, covered = {}, None
    for gene in genes:
        exps = sorted(runs.get(gene, {}),
                      key=lambda e: -len(set(runs[gene][e]) & set(structures)))
        if len(exps) < 2:
            continue
        chosen[gene] = exps[:2]
        for e in exps[:2]:
            keys = set(runs[gene][e]) & set(structures)
            covered = keys if covered is None else (covered & keys)
    if len(chosen) < 2 or covered is None or len(covered) < 100:
        return np.nan, list(chosen)
    common = sorted(covered)
    a_cols, b_cols, used = [], [], []
    for gene, (ea, eb) in chosen.items():
        a_cols.append(rankdata([runs[gene][ea][st] for st in common]))
        b_cols.append(rankdata([runs[gene][eb][st] for st in common]))
        used.append(gene)
    half = spearmanr(np.mean(a_cols, axis=0), np.mean(b_cols, axis=0)).statistic
    # each half-composite is built from half the experiments, so correct it up
    return spearman_brown(float(half)), used


def residual(y, xs):
    a = np.column_stack(list(xs) + [np.ones(len(y))])
    return y - a @ np.linalg.lstsq(a, y, rcond=None)[0]


def r2(y, xs):
    return float(1 - residual(y, xs).var() / y.var())


def splits():
    """Every way of cutting ten animals into two fives, each split once."""
    out = []
    for pick in itertools.combinations(range(len(ADULTS)), HALF):
        if 0 in pick:
            out.append((list(pick), [i for i in range(len(ADULTS)) if i not in pick]))
    return out


def spearman_brown(r):
    return 2 * r / (1 + r) if r > -1 else np.nan


def main():
    os.makedirs(OUT, exist_ok=True)
    nano = nano_per_mouse()
    expr, role = gene_profiles()

    shared = set.intersection(*[set(nano[m]) for m in ADULTS])
    shared &= set.intersection(*[set(expr[g]) for g in SUBUNITS + MARKERS])
    print('autofluorescence per structure, from the same brains:')
    auto = auto_per_mouse(shared)
    shared &= set.intersection(*[set(auto[m]) for m in ADULTS])
    structures = sorted(shared)
    ctrl = sorted(g for g in expr if role[g] == 'control_psd'
                  and all(s in expr[g] for s in structures))
    print(f'\n{len(structures)} structures in every mouse and every covariate; '
          f'{len(ctrl)} control genes usable for the PC')

    # ---------------------------------------------------------------- covariates
    abundance = composite(SUBUNITS, expr, structures)
    markers = composite(MARKERS, expr, structures)
    psd_pc1, frac = first_pc(ctrl, expr, structures)
    autofluo = rankdata([float(np.mean([auto[m][s] for m in ADULTS])) for s in structures])
    y_full = rankdata([float(np.mean([nano[m][s] for m in ADULTS])) for s in structures])
    print(f'  psd_pc1 carries {frac:.1%} of the variance of those {len(ctrl)} genes')
    print(f'\n  each covariate on its own, against the map:')
    for name, v in (('abundance (Gria1-4)', abundance), ('synaptic markers', markers),
                    ('psd_pc1', psd_pc1), ('autofluorescence', autofluo)):
        print(f'    {name:22s} rho {spearmanr(y_full, v).statistic:+.3f}')

    MODELS = [('abundance only', [abundance]),
              ('markers only', [markers]),
              ('psd_pc1 only', [psd_pc1]),
              ('autofluorescence only', [autofluo]),
              ('abundance + markers', [abundance, markers]),
              ('abundance + markers + psd_pc1', [abundance, markers, psd_pc1]),
              ('everything', [abundance, markers, psd_pc1, autofluo])]

    # ------------------------------------------------------------------ ceiling
    pairs = splits()
    raw = []
    for a, b in pairs:
        ya = rankdata([float(np.mean([nano[ADULTS[i]][s] for i in a])) for s in structures])
        yb = rankdata([float(np.mean([nano[ADULTS[i]][s] for i in b])) for s in structures])
        raw.append(spearmanr(ya, yb).statistic)
    half = float(np.mean(raw))
    ceiling = spearman_brown(half)
    print(f'\nreliability of the nano {READING} region map over {len(pairs)} splits:')
    print(f'  half-cohort against half-cohort  rho = {half:.3f}')
    print(f'  Spearman-Brown, full ten animals rho = {ceiling:.3f}   '
          f'(so at most {ceiling ** 2:.1%} of variance is explainable)')

    # ------------------------------------------- what each model explains, and what survives
    print(f'\n{"model":32s} {"R2":>7s} {"of ceiling":>11s} {"residual split-half":>20s}')
    rows, per_split = [], {}
    for label, xs in MODELS:
        explained = r2(y_full, xs)
        res_rho = []
        for a, b in pairs:
            ya = rankdata([float(np.mean([nano[ADULTS[i]][s] for i in a])) for s in structures])
            yb = rankdata([float(np.mean([nano[ADULTS[i]][s] for i in b])) for s in structures])
            res_rho.append(spearmanr(residual(ya, xs), residual(yb, xs)).statistic)
        rr = float(np.mean(res_rho))
        per_split[label] = list(res_rho)
        rows.append(dict(model=label, r2=explained,
                         share_of_ceiling=explained / max(ceiling ** 2, 1e-9),
                         residual_half=rr, residual_full=spearman_brown(rr)))
        print(f'  {label:30s} {explained:7.3f} {explained / ceiling ** 2:10.1%} '
              f'{rr:10.3f} -> {spearman_brown(rr):.3f}')

    # ---------------------------------- could better covariates close the gap?
    print(f'\nis the residual just noise in the covariates?')
    for name, genes, xs in (('abundance', SUBUNITS, [abundance]),
                            ('markers', MARKERS, [markers])):
        rel_cov, used = covariate_reliability(genes, structures)
        if np.isnan(rel_cov):
            continue
        obs = spearmanr(y_full, xs[0]).statistic
        true = obs / math.sqrt(max(ceiling * rel_cov, 1e-9))
        print(f'  {name:10s} reliability {rel_cov:.3f} (from {len(used)} of '
              f'{len(genes)} genes measured twice)')
        print(f'  {"":10s} rho {obs:+.3f} observed, {true:+.3f} corrected for '
              f'noise in both -> R2 {min(true ** 2, 1.0):.3f} against {r2(y_full, xs):.3f}')

    path = os.path.join(OUT, 'variance_partition.csv')
    with open(path, 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows({k: (f'{v:.4f}' if isinstance(v, float) else v) for k, v in r.items()}
                    for r in rows)
    print(f'\n-> {path}')

    # --------------------------------------------- where the leftover actually is
    full_model = MODELS[-1][1]
    res = residual(y_full, full_model)
    order = np.argsort(res)
    with open(os.path.join(OUT, 'residual_by_structure.csv'), 'w', newline='',
              encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(['structure', 'nano_rank', 'predicted_rank', 'residual'])
        pred = y_full - res
        for i in np.argsort(-res):
            w.writerow([structures[i], f'{y_full[i]:.1f}', f'{pred[i]:.1f}',
                        f'{res[i]:.2f}'])
    print('\nmore surface GluA1 than abundance, density and tissue predict:')
    for i in order[::-1][:10]:
        print(f'  {structures[i][:52]:54s} {res[i]:+6.1f} ranks')
    print('less:')
    for i in order[:6]:
        print(f'  {structures[i][:52]:54s} {res[i]:+6.1f} ranks')

    # is any single gene in the panel the thing we left out?
    best = []
    for gene, prof in expr.items():
        if all(s in prof for s in structures):
            rho = spearmanr(res, [prof[s] for s in structures]).statistic
            best.append((abs(rho), rho, gene))
    best.sort(reverse=True)
    print('\nthe residual against every gene in the panel -- is it just something we omitted?')
    for _, rho, gene in best[:6]:
        print(f'  {gene:10s} rho {rho:+.3f}  ({role[gene]})')

    figure(rows, ceiling, half, raw, res, structures, y_full,
           per_split[MODELS[-1][0]])


def figure(rows, ceiling, half, raw, res, structures, y_full, res_splits):
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.6),
                             gridspec_kw=dict(width_ratios=[1.25, 1.0, 1.0]))

    ax = axes[0]
    labels = [r['model'] for r in rows]
    vals = [r['r2'] for r in rows]
    y = np.arange(len(rows))
    ax.barh(y, vals, color='0.65', edgecolor='0.25', linewidth=0.5)
    ax.axvline(ceiling ** 2, color='#c0392b', lw=1.8)
    ax.annotate(f'ceiling {ceiling ** 2:.2f}\n(the map\'s own reliability)',
                (ceiling ** 2, len(rows) - 0.4), color='#c0392b', fontsize=7.5,
                ha='right', va='top', xytext=(-5, 0), textcoords='offset points')
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=7.5)
    ax.set_xlabel('variance of the nano map explained (R2, ranks)', fontsize=8)
    ax.set_title('1. neither explanation fills the map', fontsize=9)

    ax = axes[1]
    bins = np.linspace(min(min(raw), min(res_splits)) - 0.005,
                       max(max(raw), max(res_splits)) + 0.005, 30)
    ax.hist(raw, bins=bins, color='0.6', edgecolor='0.3', linewidth=0.3,
            alpha=0.85, label='the map itself')
    ax.hist(res_splits, bins=bins, color='#c0392b', edgecolor='0.3', linewidth=0.3,
            alpha=0.7, label='what is left of it')
    ax.set_xlabel('half-cohort against half-cohort (Spearman)', fontsize=8)
    ax.set_ylabel(f'splits of ten animals ({len(raw)})', fontsize=8)
    ax.legend(fontsize=7.5, frameon=False, loc='upper left')
    ax.set_title(f'2. the leftover is reproducible, not noise\n'
                 f'map {half:.3f}, residual {np.mean(res_splits):.3f}', fontsize=9)

    ax = axes[2]
    order = np.argsort(-res)
    show = list(order[:8]) + list(order[-6:])
    pos = np.arange(len(show))
    ax.barh(pos, [res[i] for i in show],
            color=['#c0392b' if res[i] > 0 else '#1f3b73' for i in show],
            edgecolor='0.25', linewidth=0.4)
    ax.set_yticks(pos)
    ax.set_yticklabels([structures[i][:34] for i in show], fontsize=6.6)
    ax.invert_yaxis()
    ax.axvline(0, color='0.3', lw=0.7)
    ax.set_xlabel('nano rank minus predicted rank', fontsize=8)
    ax.set_title('3. where the leftover lives', fontsize=9)

    for ax in axes:
        ax.tick_params(labelsize=7)
        for side in ('top', 'right'):
            ax.spines[side].set_visible(False)
    fig.suptitle('Is the nano map just receptor abundance, or just synaptic density?  '
                 'Take both out and see whether what is left survives splitting the cohort.',
                 fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    path = os.path.join(OUT, 'beyond_density.png')
    fig.savefig(path, dpi=200); plt.close(fig)
    print(f'\n{path}')


if __name__ == '__main__':
    main()
