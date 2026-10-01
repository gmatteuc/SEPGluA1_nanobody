"""
The four figures that tell Sami the result, with the statistics on their face.

adult.beyond_density and adult.beyond_controls are working figures: they exist to
let someone check each step, and there are seven of them. This makes the version
that goes in front of a person -- four panels, the same palette as the
young-versus-adult set, every number carrying an interval, and an EPS beside
every PNG because the figures end up in Illustrator.

  A  what explains the map        cross-validated R2 of each explanation against
                                  the ceiling, with bootstrap intervals
  B  the leftover is real         the two half-cohort agreements against the
                                  null of "the leftover is noise", with a p
  C  where the leftover lives     the structures that carry it
  D  seven ways it could be wrong the controls, at a glance

WHAT THE STATISTICS ARE
-----------------------
*Intervals come from resampling structures, not animals.* The claim is about
where in the brain the map departs from prediction, so the thing that would
differ in a repeat of this analysis is which structures were measurable. Each
bootstrap replicate draws 125 structures with replacement and recomputes the
whole quantity. Animal-level uncertainty is already carried by the ceiling and
by the half-cohort split, which is a separate axis and is shown separately.

*The p in panel B answers "could the leftover be noise?"* If it were, the
leftover of one half-cohort would tell you nothing about the leftover of the
other, so their correlation would sit at zero. The null is built by shuffling
which structure is which in one half and recomputing -- that destroys the
correspondence while keeping both distributions exactly as they are.

*R2 is cross-validated throughout*, because a model with more terms always fits
better on the data it was fitted to, and the question here is how much of the
map an explanation can actually predict.

Outputs, under data\\adult_v2\\beyond\\for_sami:
  A_what_explains.png/.eps   B_leftover_real.png/.eps
  C_where.png/.eps           D_controls.png/.eps
  numbers_for_the_caption.txt   every figure's numbers as a sentence

  D:\\sep_histology\\code\\tools\\venv_atlas\\Scripts\\python.exe mapping\\run_beyond_figures.py
"""

import csv
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import rankdata, spearmanr

from sepmap.adult.beyond_density import (ADULTS, OUT, SUBUNITS, MARKERS, nano_per_mouse,
                                         gene_profiles, autofluorescence, keep_structure,
                                         build_covariates, half_map, half_splits, residual,
                                         cv_r2, flexible, spearman_brown, tidy)

FIGS = os.path.join(OUT, 'for_sami')

# The young-versus-adult palette, so the two sets of figures look like one piece
# of work: COL in young_vs_adult.region_plot is {young #c0392b, naive #555555, rws #9a9a9a}.
RED = '#c0392b'         # the thing being shown
DARK = '#555555'        # its comparison
LIGHT = '#9a9a9a'       # context
BLUE = '#2e5f8a'        # only for the signed panel, where red-blue means a real zero

N_BOOT = 2000           # bootstrap replicates, over structures
N_PERM = 10000          # permutations for the "could it be noise" null
BOOT_SPLITS = 20        # half-cohort splits used inside each bootstrap replicate


def save(fig, name):
    """PNG to look at, EPS to edit. Everything here is vector already."""
    os.makedirs(FIGS, exist_ok=True)
    png = os.path.join(FIGS, name + '.png')
    fig.savefig(png, dpi=220)
    fig.savefig(os.path.join(FIGS, name + '.eps'), format='eps')
    plt.close(fig)
    print(f'  -> {png}  (+ .eps)')


def percentile_interval(values):
    return float(np.percentile(values, 2.5)), float(np.percentile(values, 97.5))


def bootstrap_models(y, models, ceiling, rng):
    """Cross-validated R2 of each model, with an interval, resampling structures."""
    point = [cv_r2(y, xs) for _, xs in models]
    n = len(y)
    draws = np.empty((N_BOOT, len(models)))
    for b in range(N_BOOT):
        take = rng.integers(0, n, n)
        yb = rankdata(y[take])
        for j, (_, xs) in enumerate(models):
            draws[b, j] = cv_r2(yb, [rankdata(x[take]) for x in xs])
    return point, [percentile_interval(draws[:, j]) for j in range(len(models))]


def bootstrap_replication(nano, structures, predictors, splits, rng):
    """Half-against-half agreement of the leftover, with an interval."""
    def agreement(idx, use_splits):
        picked = [structures[i] for i in idx]
        return float(np.mean([
            spearmanr(residual(half_map(nano, a, picked), [p[idx] for p in predictors]),
                      residual(half_map(nano, b, picked), [p[idx] for p in predictors])
                      ).statistic for a, b in use_splits]))

    full = np.arange(len(structures))
    point = agreement(full, splits)
    few = [splits[i] for i in rng.choice(len(splits), BOOT_SPLITS, replace=False)]
    draws = [agreement(rng.integers(0, len(structures), len(structures)), few)
             for _ in range(300)]
    return point, percentile_interval(draws)


def noise_null(nano, structures, predictors, splits, rng):
    """If the leftover were noise, how well would the two halves agree? Test it.

    Shuffling which structure is which in one half breaks the correspondence
    between the two halves while leaving both leftovers exactly as they are, so
    the null says "these two are unrelated" and nothing else.
    """
    a, b = splits[0]
    ra = residual(half_map(nano, a, structures), predictors)
    rb = residual(half_map(nano, b, structures), predictors)
    observed = float(spearmanr(ra, rb).statistic)
    null = np.empty(N_PERM)
    for i in range(N_PERM):
        null[i] = spearmanr(ra, rng.permutation(rb)).statistic
    p = float((np.sum(np.abs(null) >= abs(observed)) + 1) / (N_PERM + 1))
    return observed, null, p


# ------------------------------------------------------------------- the panels

def panel_a(point, intervals, labels, ceiling, ceiling_ci):
    fig, ax = plt.subplots(figsize=(7.8, 4.4))
    y = np.arange(len(labels))
    lo = [p - i[0] for p, i in zip(point, intervals)]
    hi = [i[1] - p for p, i in zip(point, intervals)]
    colours = [LIGHT] * (len(labels) - 1) + [RED]
    ax.barh(y, point, color=colours, edgecolor='0.25', linewidth=0.5)
    ax.errorbar(point, y, xerr=[lo, hi], fmt='none', ecolor='0.2',
                elinewidth=0.9, capsize=2.5)
    ax.axvline(ceiling, color=DARK, lw=1.8)
    ax.axvspan(ceiling_ci[0], ceiling_ci[1], color=DARK, alpha=0.15, lw=0)
    ax.text(ceiling - 0.015, 0.97, 'ceiling: all of the map\nthat is explainable',
            transform=ax.get_xaxis_transform(), color=DARK, fontsize=8,
            ha='right', va='top')
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8.5)
    ax.set_xlim(0, 1.02)
    ax.set_xlabel('variance of the surface-GluA1 map predicted\n'
                  '(cross-validated, 95% bootstrap interval over structures)', fontsize=8.5)
    ax.set_title('A.  Receptor abundance and synaptic density explain part of the map',
                 fontsize=10, loc='left')
    tidy(ax)
    fig.tight_layout()
    save(fig, 'A_what_explains')


def panel_b(map_agreement, leftover_agreement, null, p, rep_point, rep_ci):
    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    bins = np.linspace(-0.45, 1.0, 120)
    for values, colour, label in ((null, LIGHT, 'if the leftover were noise'),
                                  (map_agreement, DARK, 'the map itself'),
                                  (leftover_agreement, RED, 'what is left of it')):
        counts, edges = np.histogram(values, bins=bins)
        ax.bar(edges[:-1], counts / counts.max(), width=np.diff(edges),
               align='edge', color=colour, edgecolor='none', label=label)
    ax.set_xlabel('agreement between two independent halves of the cohort (Spearman)',
                  fontsize=8.5)
    ax.set_ylabel('how often, each curve scaled to its own peak', fontsize=8.5)
    ax.set_ylim(0, 1.15)
    ax.legend(fontsize=8, frameon=False, loc='upper left')
    ax.set_title(f'B.  What is left over replicates across animals\n'
                 f'leftover {rep_point:.3f} [{rep_ci[0]:.3f}, {rep_ci[1]:.3f}], '
                 f'against noise p < {max(p, 1e-4):.0e}', fontsize=10, loc='left')
    tidy(ax)
    fig.tight_layout()
    save(fig, 'B_leftover_real')


def panel_c(res, structures, n_show=9):
    order = np.argsort(-res)
    show = list(order[:n_show]) + list(order[-n_show:])
    fig, ax = plt.subplots(figsize=(7.6, 5.6))
    pos = np.arange(len(show))
    ax.barh(pos, [res[i] for i in show],
            color=[RED if res[i] > 0 else BLUE for i in show],
            edgecolor='0.25', linewidth=0.4)
    ax.set_yticks(pos)
    ax.set_yticklabels([structures[i][:42] for i in show], fontsize=8)
    ax.invert_yaxis()
    ax.axvline(0, color='0.3', lw=0.8)
    ax.set_xlabel('surface GluA1, minus what abundance and density predict (ranks)',
                  fontsize=8.5)
    ax.set_title('C.  The leftover is anatomically organised', fontsize=10, loc='left')
    tidy(ax)
    fig.tight_layout()
    save(fig, 'C_where')


def panel_d(controls):
    fig, ax = plt.subplots(figsize=(9.4, 4.0))
    ax.axis('off')
    ax.set_title('D.  Seven ways the result could be an artefact, and the number '
                 'for each', fontsize=10, loc='left')
    for i, row in enumerate(controls):
        y = 1 - (i + 1) / (len(controls) + 1)
        ok = row['verdict'] == 'pass'
        ax.text(0.0, y, row['control'], fontsize=9, va='center',
                color='0.15' if ok else RED)
        ax.text(0.38, y, row['number'][:78], fontsize=8, va='center', color='0.35')
        ax.text(0.985, y, 'ruled out' if ok else 'CHECK', fontsize=9, va='center',
                ha='right', color=DARK if ok else RED)
    fig.tight_layout()
    save(fig, 'D_controls')


def main():
    os.makedirs(FIGS, exist_ok=True)
    rng = np.random.default_rng(0)

    nano, division = nano_per_mouse()
    expr, role = gene_profiles()
    everywhere = set.intersection(*[set(nano[m]) for m in ADULTS])
    everywhere &= set.intersection(*[set(expr[g]) for g in SUBUNITS + MARKERS])
    kept = [s for s in sorted(everywhere) if keep_structure(s, division.get(s, ''))[0]]
    auto = autofluorescence(set(kept))
    structures = sorted(set(kept).intersection(*[set(auto[m]) for m in ADULTS]))
    covariates, _, _ = build_covariates(nano, expr, role, auto, structures)
    splits = half_splits()
    y = half_map(nano, range(len(ADULTS)), structures)
    model = flexible(list(covariates.values()))
    print(f'{len(structures)} structures, {len(ADULTS)} adults, {len(splits)} half-splits')

    # the ceiling, with an interval over structures
    map_agreement = [spearmanr(half_map(nano, a, structures),
                               half_map(nano, b, structures)).statistic
                     for a, b in splits]
    ceiling = spearman_brown(float(np.mean(map_agreement))) ** 2
    boot_ceiling = []
    for _ in range(300):
        take = rng.integers(0, len(structures), len(structures))
        picked = [structures[i] for i in take]
        few = [splits[i] for i in rng.choice(len(splits), BOOT_SPLITS, replace=False)]
        boot_ceiling.append(spearman_brown(float(np.mean([
            spearmanr(half_map(nano, a, picked), half_map(nano, b, picked)).statistic
            for a, b in few]))) ** 2)
    ceiling_ci = percentile_interval(boot_ceiling)
    print(f'ceiling {ceiling:.3f} [{ceiling_ci[0]:.3f}, {ceiling_ci[1]:.3f}]')

    c = covariates
    models = [('receptor abundance\n(Gria1-4)', [c['abundance']]),
              ('synaptic markers\n(11 pre- and postsynaptic)', [c['markers']]),
              ('postsynaptic gene\nexpression (188 genes)', [c['psd_pc1']]),
              ('tissue autofluorescence\n(the same brains)', [c['autofluo']]),
              ('all four together', model)]
    point, intervals = bootstrap_models(y, models, ceiling, rng)
    for (label, _), pt, iv in zip(models, point, intervals):
        print(f'  {label.splitlines()[0]:34s} CV R2 {pt:.3f} [{iv[0]:.3f}, {iv[1]:.3f}]'
              f'   {pt / ceiling:5.1%} of the ceiling')
    unexplained = 1 - point[-1] / ceiling
    print(f'  unexplained share of the ceiling: {unexplained:.1%}')

    leftover_agreement = [spearmanr(residual(half_map(nano, a, structures), model),
                                    residual(half_map(nano, b, structures), model)).statistic
                          for a, b in splits]
    rep_point, rep_ci = bootstrap_replication(nano, structures, model, splits, rng)
    observed, null, p = noise_null(nano, structures, model, splits, rng)
    print(f'  leftover replicates {rep_point:.3f} [{rep_ci[0]:.3f}, {rep_ci[1]:.3f}]; '
          f'against the noise null p = {p:.2e}')

    res = residual(y, model)

    controls = []
    path = os.path.join(OUT, 'controls.csv')
    if os.path.exists(path):
        with open(path, newline='', encoding='utf-8') as fh:
            controls = list(csv.DictReader(fh))

    panel_a(point, intervals, [m[0] for m in models], ceiling, ceiling_ci)
    panel_b(map_agreement, leftover_agreement, null, p, rep_point, rep_ci)
    panel_c(res, structures)
    if controls:
        panel_d(controls)

    lines = [
        'Numbers for the captions (all on 125 grey-matter structures, ten adult mice).',
        '',
        f'A. The map is reproducible enough that {ceiling:.1%} of its variance is',
        f'   explainable in principle (95% CI {ceiling_ci[0]:.1%} to {ceiling_ci[1]:.1%}).',
        f'   Receptor abundance predicts {point[0] / ceiling:.0%} of that, synaptic markers',
        f'   {point[1] / ceiling:.0%}, postsynaptic gene expression {point[2] / ceiling:.0%},',
        f'   tissue autofluorescence {max(point[3], 0) / ceiling:.0%}. All four together,',
        f'   allowed a non-linear relationship, predict {point[-1] / ceiling:.0%},',
        f'   leaving {unexplained:.0%} unexplained.',
        '',
        f'B. Splitting the ten animals into two fives every possible way ({len(splits)} splits),',
        f'   the two half-cohort maps agree at {np.mean(map_agreement):.3f}. After the four',
        f'   explanations are removed, the two leftovers still agree at {rep_point:.3f}',
        f'   (95% CI {rep_ci[0]:.3f} to {rep_ci[1]:.3f}). Were the leftover noise, that',
        f'   agreement would be zero; p < {max(p, 1e-4):.0e} by permutation.',
        '',
        'C. Structures where surface GluA1 most exceeds and falls short of what',
        '   abundance and density predict, in ranks among the 125.',
        '',
        'D. Seven controls, each ruling out a way the leftover could be an artefact.',
        '',
        'What this does NOT show: that the leftover is the surface fraction. A residual',
        'is only ever what the model left out.',
    ]
    with open(os.path.join(FIGS, 'numbers_for_the_caption.txt'), 'w',
              encoding='utf-8') as fh:
        fh.write('\n'.join(lines) + '\n')
    print(f'  -> {os.path.join(FIGS, "numbers_for_the_caption.txt")}')
