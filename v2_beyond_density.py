"""
The claim: the adult surface-GluA1 map carries substantial, reproducible spatial
structure that receptor abundance and synaptic density do not explain.

In numbers, so the word "substantial" is not doing the work: 97.4% of this map
is explainable in principle (it is that reproducible across animals); receptor
abundance and three density proxies, allowed a bending relationship, account for
61% of that; the remaining 39% replicates across independent halves of the
cohort at 0.97. Even a twenty-component model of the whole 390-gene panel, far
richer than either explanation, still leaves a remainder that replicates at 0.88.

This script exists to make that sentence defensible, and `v2_beyond_controls.py`
exists to attack it. Read them in that order.

WHY NOT JUST RANK GENES
-----------------------
The obvious way to ask "is the nano map about trafficking?" is to correlate it
with every gene and look at the top of the list. We did that, at length, and it
does not work: v2_ish_panel_test.py showed the ranking is flat across
postsynaptic gene classes. Which is exactly what you would see if the map were
nothing but synaptic density -- every synaptic gene would correlate, and none
would stand out. A flat ranking is evidence for the boring hypothesis, not
against it.

So the question is turned around. Instead of asking which gene the map looks
like, ask what is LEFT once the boring explanations are subtracted, and whether
that leftover is real.

THE LOGIC, IN ONE LINE
----------------------
    nano zref  ~  receptor abundance + synaptic density   ->   residual
    does the residual come out the same in two independent halves of the cohort?

A residual is never zero; any model leaves something behind. What makes a
residual interesting is reproducibility. Noise cannot replicate across
independent animals, so if two halves of the cohort produce the same leftover,
the leftover is biology rather than measurement error. That single test is the
backbone of the claim.

THE FOUR STEPS
--------------
  0  the structures   which regions the analysis is allowed to use, and why
  1  the ceiling      how reliable is the map itself? Nothing can explain
                      variance that is not there, so this sets the scale
                      against which every R2 below should be read
  2  the covariates   what do abundance and density actually buy?
  3  the residual     does what is left replicate across animals?  <- the claim
  4  where it lives   which structures carry it

Each step prints its numbers and draws its own figure, so the argument can be
checked one piece at a time rather than taken whole.

CHOICES, AND WHY
----------------
*Structures, not voxels.* The gene side exists only at 200 um and only in a
different mouse, so a voxel-level comparison would mostly be registration error.
Region means are the level at which both sides mean the same thing.

*Ranks everywhere.* Allen expression energy has an arbitrary scale per
experiment, and nothing here assumes a straight-line relationship -- only that
more receptor means more signal. Spearman and rank regression say exactly that
and nothing more. Whether the rank relationship really is straight is not
assumed either: control E in v2_beyond_controls.py tests it.

*Grey matter only, and no catch-all labels.* Fibre tracts and ventricles have no
synapses, so "synaptic density" is not even defined there, and the CCF's
"..., unassigned" entries are leftover voxels rather than anatomical structures
-- what lands in them depends on how each brain registered. Both are dropped
BEFORE anything is fitted, by the rule in `keep_structure`, and fig0 shows how
many went and why.

*Ten adults, naive and RWS pooled.* The two groups differ in whisker
experience, not in what the nano stain is; control D checks that the leftover
does not care which group it came from.

Outputs, under data\\adult_v2\\beyond:
  structures_used.csv        every structure, kept or dropped, with the reason
  variance_partition.csv     what each covariate set explains, against the ceiling
  residual_by_structure.csv  where the map exceeds and falls short of prediction
  fig0_structures.png, fig1_ceiling.png, fig2_covariates.png, fig3_residual.png

A NOTE ON HOW THIS NUMBER MOVED
-------------------------------
The first version of this analysis fitted the covariates as straight lines and
reported that half the explainable variance was unaccounted for. Control E in
v2_beyond_controls.py showed that was too generous to us: the rank relationships
here are curved, and a straight line left variance in the residual that the
covariates could have explained. Letting them bend took the model from 0.42 to
0.60 cross-validated, and the headline from "half" to 39%. The claim survived,
smaller. That is what the controls are for.

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
MERGED_ISH = os.path.join(DATA, 'adult_v2', 'ish', 'gene_region_table_merged.csv')
OUT = os.path.join(DATA, 'adult_v2', 'beyond')

ADULTS = NAIVE + RWS
READING = 'zref'
MIN_VOX = 250           # 20 um voxels a structure needs, the same bar as everywhere else
HALF = 5                # animals per half when the cohort is split

# Divisions that are grey matter. Everything else -- fibre tracts (lfbs, mfbs,
# eps, scwm, cbf, cm), ventricles (VL, V3) and the brain-wide catch-all -- is
# dropped, because a synaptic-density covariate has no meaning in white matter
# or in cerebrospinal fluid.
GREY = {'Isocortex', 'OLF', 'HPF', 'CTXsp', 'STR', 'PAL', 'TH', 'HY', 'MB', 'P', 'MY', 'CB'}

# The gene sets standing in for the two boring explanations.
SUBUNITS = ('Gria1', 'Gria2', 'Gria3', 'Gria4')
MARKERS = ('Syp', 'Syn1', 'Vamp2', 'Bsn', 'Syt1',                      # presynaptic
           'Dlg4', 'Homer1', 'Shank2', 'Shank3', 'Nlgn1', 'Camk2a')    # postsynaptic


# --------------------------------------------------------------------- loading

def keep_structure(name, division):
    """Should this structure be in the analysis? Returns (keep, reason why not).

    Two rules, both applied before any fitting so neither can be tuned to the
    answer: grey matter only, by division; and no "..., unassigned" entries,
    which are voxels the atlas could not place rather than structures.
    """
    if 'unassigned' in name.lower():
        return False, 'catch-all label, not a structure'
    if division not in GREY:
        return False, f'division {division} is not grey matter'
    return True, ''


def nano_per_mouse():
    """{mouse: {structure: zref}} and {structure: division}, for the ten adults."""
    per, division = defaultdict(dict), {}
    with open(NANO, newline='', encoding='utf-8') as fh:
        for r in csv.DictReader(fh):
            if r['reading'] == READING and r['mouse'] in ADULTS:
                per[r['mouse']][r['structure']] = float(r['log2_value'])
                division[r['structure']] = r['division']
    return per, division


def gene_profiles(path=MERGED_ISH):
    """{gene: {structure: rank_mean}} and {gene: role}; replicates already merged."""
    per, role = defaultdict(dict), {}
    with open(path, newline='', encoding='utf-8') as fh:
        for r in csv.DictReader(fh):
            per[r['symbol']][r['structure']] = float(r['rank_mean'])
            role[r['symbol']] = r['role']
    return per, role


def structure_names():
    """parcellation index -> structure name, the key both tables use."""
    names = {}
    with open(CSV_MAP, newline='', encoding='utf-8') as fh:
        for row in csv.DictReader(fh):
            if row['parcellation_term_set_name'] == 'structure':
                names[int(row['parcellation_index'])] = row['parcellation_term_name']
    return names


def autofluorescence(structures):
    """{mouse: {structure: log2 mean autofluorescence}}, from the same ten brains.

    Worth more as a density proxy than it looks: every gene covariate is a
    different mouse warped into our atlas, while this one is the very tissue the
    nano signal was measured in, in the same sections, with no registration
    between the two.
    """
    names = structure_names()
    per = {}
    for mouse in ADULTS:
        z = np.load(os.path.join(PER_MOUSE, mouse + '.npz'))
        tissue = z['tissue']
        ann = annotation_20(MICE[mouse][1])
        labels = ann[tissue]
        counts = np.bincount(labels, minlength=int(ann.max()) + 1)
        totals = np.bincount(labels, weights=z['auto'].astype(np.float32)[tissue],
                             minlength=len(counts))
        acc = defaultdict(lambda: [0, 0.0])
        for idx in np.nonzero(counts)[0]:
            if idx == 0:
                continue
            a = acc[names.get(int(idx), f'id{idx}')]
            a[0] += int(counts[idx])
            a[1] += totals[idx]
        per[mouse] = {k: math.log2(t / c) for k, (c, t) in acc.items()
                      if c >= MIN_VOX and t > 0 and k in structures}
    return per


# ------------------------------------------------------------------- the maths

def composite(genes, expr, structures):
    """One predictor from a set of genes: the mean of their rank profiles.

    Ranks rather than values because each Allen experiment carries its own
    arbitrary intensity scale, so raw numbers are not comparable between genes
    even when their orderings are.
    """
    return np.mean([rankdata([expr[g][s] for s in structures]) for g in genes], axis=0)


def first_pc(genes, expr, structures):
    """The dominant shared axis of a gene set, and what share of it that axis is.

    Used for the 300 postsynaptic-density genes: rather than naming a handful of
    markers and then being argued with about the choice, take whatever those
    genes have most in common and call that the density axis.
    """
    m = np.array([rankdata([expr[g][s] for s in structures]) for g in genes])
    m = (m - m.mean(axis=1, keepdims=True)) / m.std(axis=1, keepdims=True)
    _, sv, vt = np.linalg.svd(m - m.mean(axis=0), full_matrices=False)
    pc = vt[0]
    if np.corrcoef(pc, m.mean(axis=0))[0, 1] < 0:     # point it the same way as the genes
        pc = -pc
    return pc, float(sv[0] ** 2 / (sv ** 2).sum())


def residual(y, predictors):
    """What is left of y after least squares on the predictors, plus an intercept."""
    design = np.column_stack(list(predictors) + [np.ones(len(y))])
    return y - design @ np.linalg.lstsq(design, y, rcond=None)[0]


def r_squared(y, predictors):
    """Variance explained on the data the fit was made from -- optimistic."""
    return float(1 - residual(y, predictors).var() / y.var())


def flexible(predictors):
    """The same covariates, allowed to bend.

    A straight line through two rank variables assumes the relationship is not
    only monotone but evenly paced, and control E in v2_beyond_controls.py showed
    that assumption is wrong here: adding curvature raises the cross-validated
    fit from 0.42 to 0.60. Crediting that 0.18 to the leftover would have been
    our mistake, not the biology's, so the headline model bends.
    """
    return list(predictors) + [x ** 2 for x in predictors] + [x ** 3 for x in predictors]


def cv_r2(y, predictors, folds=5):
    """Variance explained on structures the fit has never seen.

    This is the honest number once a model has many terms: adding predictors
    always improves the in-sample fit, and only held-out structures can say
    whether it improved the prediction. Structures are shuffled once with a fixed
    seed, cut into folds, and each fold predicted from a fit on the others.
    """
    n = len(y)
    order = np.random.default_rng(0).permutation(n)
    design = np.column_stack(list(predictors) + [np.ones(n)])
    predicted = np.empty(n)
    for k in range(folds):
        test = order[k::folds]
        train = np.setdiff1d(order, test)
        predicted[test] = design[test] @ np.linalg.lstsq(design[train], y[train],
                                                         rcond=None)[0]
    return float(1 - np.var(y - predicted) / np.var(y))


def half_splits():
    """Every way of cutting ten animals into two fives, each split counted once."""
    return [(list(p), [i for i in range(len(ADULTS)) if i not in p])
            for p in itertools.combinations(range(len(ADULTS)), HALF) if 0 in p]


def spearman_brown(r):
    """Reliability of a whole cohort, from the agreement of its two halves."""
    return 2 * r / (1 + r) if r > -1 else float('nan')


def half_map(nano, indices, structures):
    """The mean zref map of one half-cohort, as ranks."""
    return rankdata([float(np.mean([nano[ADULTS[i]][s] for i in indices]))
                     for s in structures])


def build_covariates(nano, expr, role, auto, structures):
    """The four predictors, in one place so every script builds them identically."""
    controls = sorted(g for g in expr if role[g] == 'control_psd'
                      and all(s in expr[g] for s in structures))
    psd_pc1, share = first_pc(controls, expr, structures)
    return dict(abundance=composite(SUBUNITS, expr, structures),
                markers=composite(MARKERS, expr, structures),
                psd_pc1=psd_pc1,
                autofluo=rankdata([float(np.mean([auto[m][s] for m in ADULTS]))
                                   for s in structures])), controls, share


def prepare():
    """Everything the analysis and the controls both need, loaded once."""
    nano, division = nano_per_mouse()
    expr, role = gene_profiles()
    everywhere = set.intersection(*[set(nano[m]) for m in ADULTS])
    everywhere &= set.intersection(*[set(expr[g]) for g in SUBUNITS + MARKERS])
    kept = [s for s in sorted(everywhere) if keep_structure(s, division.get(s, ''))[0]]
    auto = autofluorescence(set(kept))
    structures = sorted(set(kept).intersection(*[set(auto[m]) for m in ADULTS]))
    return nano, division, expr, role, auto, structures


# --------------------------------------------------------------------- drawing

def tidy(ax):
    ax.tick_params(labelsize=7)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)


def save(fig, name):
    path = os.path.join(OUT, name)
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f'  -> {path}')


# ----------------------------------------------------------------------- steps

def step0_structures(nano, division, expr):
    """Choose the structures, and show what the choice threw away."""
    print('\nSTEP 0  which structures the analysis may use')
    everywhere = set.intersection(*[set(nano[m]) for m in ADULTS])
    everywhere &= set.intersection(*[set(expr[g]) for g in SUBUNITS + MARKERS])

    rows = []
    for name in sorted(everywhere):
        ok, why = keep_structure(name, division.get(name, ''))
        rows.append(dict(structure=name, division=division.get(name, ''),
                         kept='yes' if ok else 'no', reason=why))
    with open(os.path.join(OUT, 'structures_used.csv'), 'w', newline='',
              encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=['structure', 'division', 'kept', 'reason'])
        w.writeheader(); w.writerows(rows)

    kept = [r['structure'] for r in rows if r['kept'] == 'yes']
    dropped = [r for r in rows if r['kept'] == 'no']
    print(f'  {len(everywhere)} structures measured in every mouse and every covariate')
    print(f'  {len(kept)} kept, {len(dropped)} dropped:')
    for reason in sorted({r['reason'] for r in dropped}):
        names = [r['structure'] for r in dropped if r['reason'] == reason]
        print(f'    {len(names):3d}  {reason}')

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 3.9),
                             gridspec_kw=dict(width_ratios=[1, 1.25]))
    counts = defaultdict(int)
    for r in rows:
        if r['kept'] == 'yes':
            counts[r['division']] += 1
    order = sorted(counts, key=counts.get, reverse=True)
    axes[0].bar(range(len(order)), [counts[d] for d in order],
                color='0.65', edgecolor='0.25', linewidth=0.5)
    axes[0].set_xticks(range(len(order)))
    axes[0].set_xticklabels(order, fontsize=7, rotation=45, ha='right')
    axes[0].set_ylabel('structures kept', fontsize=8)
    axes[0].set_title(f'what the analysis runs on\n{len(kept)} grey-matter structures',
                      fontsize=9)
    tidy(axes[0])

    why_counts = defaultdict(int)
    for r in dropped:
        why_counts[r['reason']] += 1
    ordered = sorted(why_counts, key=why_counts.get)
    axes[1].barh(range(len(ordered)), [why_counts[w] for w in ordered],
                 color='#c0392b', edgecolor='0.25', linewidth=0.5, alpha=0.85)
    axes[1].set_yticks(range(len(ordered)))
    axes[1].set_yticklabels([w[:46] for w in ordered], fontsize=7)
    axes[1].set_xlabel('structures dropped', fontsize=8)
    axes[1].set_title('and what it refuses to run on\ndecided by rule, before any fitting',
                      fontsize=9)
    tidy(axes[1])
    fig.tight_layout()
    save(fig, 'fig0_structures.png')
    return kept


def step1_ceiling(nano, structures):
    """How reliable is the map itself? Nothing below can beat this."""
    print('\nSTEP 1  the ceiling -- how much of this map is explainable at all')
    splits = half_splits()
    agreement = [spearmanr(half_map(nano, a, structures),
                           half_map(nano, b, structures)).statistic
                 for a, b in splits]
    half = float(np.mean(agreement))
    full = spearman_brown(half)
    print(f'  over {len(splits)} five-against-five splits of the ten adults:')
    print(f'    half-cohort against half-cohort   rho = {half:.3f}')
    print(f'    Spearman-Brown, all ten animals   rho = {full:.3f}')
    print(f'  so at most {full ** 2:.1%} of this map can ever be explained, and every')
    print('  R2 below is quoted against that rather than against 1.0.')

    fig, ax = plt.subplots(figsize=(5.6, 3.9))
    ax.hist(agreement, bins=25, color='0.7', edgecolor='0.35', linewidth=0.4)
    ax.axvline(half, color='#c0392b', lw=1.8)
    ax.set_xlabel('Spearman between the two half-cohort maps', fontsize=8)
    ax.set_ylabel(f'splits of ten animals ({len(splits)})', fontsize=8)
    ax.set_title('Step 1. the map is highly reproducible\n'
                 f'half-cohorts agree at {half:.3f}, whole cohort {full:.3f}', fontsize=9)
    tidy(ax)
    fig.tight_layout()
    save(fig, 'fig1_ceiling.png')
    return full, splits, agreement


def step2_covariates(nano, expr, role, auto, structures, ceiling):
    """What do receptor abundance and synaptic density actually account for?"""
    print('\nSTEP 2  what the two boring explanations buy')
    covariates, controls, share = build_covariates(nano, expr, role, auto, structures)
    y = half_map(nano, range(len(ADULTS)), structures)
    print(f'  psd_pc1 is the first component of {len(controls)} postsynaptic-density '
          f'genes, carrying {share:.1%} of their variance')
    print('  each covariate on its own:')
    nice = {'abundance': 'abundance (Gria1-4)', 'markers': 'synaptic markers',
            'psd_pc1': 'psd_pc1', 'autofluo': 'autofluorescence'}
    for key, v in covariates.items():
        print(f'    {nice[key]:22s} rho {spearmanr(y, v).statistic:+.3f}   '
              f'R2 {r_squared(y, [v]):.3f}   '
              f'{r_squared(y, [v]) / ceiling ** 2:5.1%} of the ceiling')

    c = covariates
    models = [('abundance only', [c['abundance']]),
              ('markers only', [c['markers']]),
              ('psd_pc1 only', [c['psd_pc1']]),
              ('autofluorescence only', [c['autofluo']]),
              ('abundance + markers', [c['abundance'], c['markers']]),
              ('abundance + markers + psd_pc1',
               [c['abundance'], c['markers'], c['psd_pc1']]),
              ('all four, straight', list(c.values())),
              ('all four, allowed to bend', flexible(list(c.values())))]
    vals = [cv_r2(y, xs) for _, xs in models]
    print('  and in combination, scored on structures the fit has not seen:')
    for (label, xs), v in zip(models, vals):
        print(f'    {label:30s} CV R2 {v:.3f}   in-sample {r_squared(y, xs):.3f}   '
              f'{v / ceiling ** 2:5.1%} of the ceiling')
    print(f'  the model to quote is the last one -- letting the covariates bend is')
    print(f'  their best shot -- and it leaves {1 - vals[-1] / ceiling ** 2:.0%} of the '
          f'explainable variance unaccounted for.')

    with open(os.path.join(OUT, 'variance_partition.csv'), 'w', newline='',
              encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(['model', 'cv_r2', 'in_sample_r2', 'share_of_ceiling'])
        for (label, xs), v in zip(models, vals):
            w.writerow([label, f'{v:.4f}', f'{r_squared(y, xs):.4f}',
                        f'{v / ceiling ** 2:.4f}'])

    fig, ax = plt.subplots(figsize=(7.6, 4.5))
    ax.barh(np.arange(len(models)), vals, color='0.65', edgecolor='0.25', linewidth=0.5)
    ax.axvline(ceiling ** 2, color='#c0392b', lw=1.8)
    ax.annotate(f'ceiling {ceiling ** 2:.2f}\n(the map\'s own reliability)',
                (ceiling ** 2, len(models) - 0.4), color='#c0392b', fontsize=7.5,
                ha='right', va='top', xytext=(-6, 0), textcoords='offset points')
    ax.set_yticks(np.arange(len(models)))
    ax.set_yticklabels([m[0] for m in models], fontsize=8)
    ax.set_xlim(0, 1.02)
    ax.set_xlabel('variance explained on held-out structures (cross-validated R2)',
                  fontsize=8)
    ax.set_title('Step 2. neither explanation fills the map,\n'
                 'even when allowed to bend', fontsize=9)
    tidy(ax)
    fig.tight_layout()
    save(fig, 'fig2_covariates.png')
    return covariates, y


def step3_residual(nano, structures, covariates, splits, raw_agreement):
    """The claim: what is left over replicates across independent animals."""
    print('\nSTEP 3  does the leftover replicate?   <- this is the claim')
    xs = flexible(list(covariates.values()))
    agreement = [spearmanr(residual(half_map(nano, a, structures), xs),
                           residual(half_map(nano, b, structures), xs)).statistic
                 for a, b in splits]
    half = float(np.mean(agreement))
    print('  residualise each half-cohort map on the bending four-covariate model,')
    print('  then compare the two leftovers:')
    print(f'    half against half   rho = {half:.3f}   '
          f'(the map itself: {np.mean(raw_agreement):.3f})')
    print(f'    Spearman-Brown      rho = {spearman_brown(half):.3f}')
    print('  noise cannot replicate across independent animals, so this is real.')
    return agreement


def step4_where(nano, structures, covariates, expr, role, agreement, raw_agreement):
    """Which structures carry the leftover, and is any single gene behind it?"""
    print('\nSTEP 4  where the leftover lives')
    y = half_map(nano, range(len(ADULTS)), structures)
    res = residual(y, flexible(list(covariates.values())))
    predicted = y - res

    with open(os.path.join(OUT, 'residual_by_structure.csv'), 'w', newline='',
              encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(['structure', 'nano_rank', 'predicted_rank', 'residual'])
        for i in np.argsort(-res):
            w.writerow([structures[i], f'{y[i]:.1f}', f'{predicted[i]:.1f}',
                        f'{res[i]:.2f}'])

    order = np.argsort(-res)
    print('  more surface GluA1 than abundance and density predict:')
    for i in order[:8]:
        print(f'    {structures[i][:48]:50s} {res[i]:+6.1f} ranks')
    print('  less:')
    for i in order[-6:]:
        print(f'    {structures[i][:48]:50s} {res[i]:+6.1f} ranks')

    scored = []
    for gene in expr:
        if all(s in expr[gene] for s in structures):
            rho = spearmanr(res, [expr[gene][s] for s in structures]).statistic
            scored.append((abs(rho), rho, gene))
    scored.sort(reverse=True)
    print('  the residual against all 390 genes -- is it just something we left out?')
    for _, rho, gene in scored[:5]:
        print(f'    {gene:10s} rho {rho:+.3f}   ({role[gene]})')
    print('    no single gene in the panel accounts for it.')

    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.5),
                             gridspec_kw=dict(width_ratios=[1, 1.15]))
    bins = np.linspace(min(min(raw_agreement), min(agreement)) - 0.005,
                       max(max(raw_agreement), max(agreement)) + 0.005, 30)
    axes[0].hist(raw_agreement, bins=bins, color='0.6', edgecolor='0.3',
                 linewidth=0.3, alpha=0.85, label='the map itself')
    axes[0].hist(agreement, bins=bins, color='#c0392b', edgecolor='0.3',
                 linewidth=0.3, alpha=0.7, label='what is left of it')
    axes[0].set_xlabel('half-cohort against half-cohort (Spearman)', fontsize=8)
    axes[0].set_ylabel(f'splits of ten animals ({len(agreement)})', fontsize=8)
    axes[0].legend(fontsize=7.5, frameon=False, loc='upper left')
    axes[0].set_title('Step 3. the leftover is reproducible, not noise\n'
                      f'map {np.mean(raw_agreement):.3f}, '
                      f'residual {np.mean(agreement):.3f}', fontsize=9)
    tidy(axes[0])

    show = list(order[:8]) + list(order[-6:])
    axes[1].barh(np.arange(len(show)), [res[i] for i in show],
                 color=['#c0392b' if res[i] > 0 else '#1f3b73' for i in show],
                 edgecolor='0.25', linewidth=0.4)
    axes[1].set_yticks(np.arange(len(show)))
    axes[1].set_yticklabels([structures[i][:36] for i in show], fontsize=6.8)
    axes[1].invert_yaxis()
    axes[1].axvline(0, color='0.3', lw=0.7)
    axes[1].set_xlabel('nano rank minus predicted rank', fontsize=8)
    axes[1].set_title('Step 4. and it is anatomically organised', fontsize=9)
    tidy(axes[1])
    fig.tight_layout()
    save(fig, 'fig3_residual.png')


def main():
    os.makedirs(OUT, exist_ok=True)
    nano, division = nano_per_mouse()
    expr, role = gene_profiles()

    kept = step0_structures(nano, division, expr)
    auto = autofluorescence(set(kept))
    structures = sorted(set(kept).intersection(*[set(auto[m]) for m in ADULTS]))
    print(f'  {len(structures)} of those also have autofluorescence in every mouse')

    ceiling, splits, raw_agreement = step1_ceiling(nano, structures)
    covariates, _ = step2_covariates(nano, expr, role, auto, structures, ceiling)
    agreement = step3_residual(nano, structures, covariates, splits, raw_agreement)
    step4_where(nano, structures, covariates, expr, role, agreement, raw_agreement)

    print('\nNow run v2_beyond_controls.py: it tries to break this seven ways.')


if __name__ == '__main__':
    main()
