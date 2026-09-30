"""
The regression itself, shown as a regression: observed, predicted, residual.

Everything else in this folder reports the regression's summary numbers -- R2,
replication, controls. None of them shows the fit. This does, three ways, because
a reader's first fair question is "show me the actual prediction":

  1  observed against predicted, one dot per structure, with the identity line.
     If the model were good the cloud would sit on that line; the spread away
     from it IS the leftover, drawn rather than summarised.
  2  residual against predicted, the standard diagnostic. A tilt or a fan here
     would mean the model is mis-specified rather than merely incomplete.
  3  the same three quantities painted back onto the brain -- observed map,
     predicted map, residual map -- because the leftover is a spatial claim and
     ought to be looked at spatially.

WHAT IS BEING REGRESSED ON WHAT
-------------------------------
    y            the ten adults' mean surface-GluA1 (zref) per structure, ranked
    predictors   four, each also one value per structure, each entered as
                 x, x^2, x^3:

      abundance  the four AMPA receptor subunit genes Gria1-4, averaged over
                 their rank profiles. The set is not ours: it is GO:0004971
                 intersected with GO:0032281, less the delta receptors Grid1
                 and Grid2 (see OVERRIDE in v2_panel_build.py).
      markers    eleven canonical synaptic markers CHOSEN BY HAND from the
                 panel -- Syp, Syn1, Vamp2, Bsn, Syt1 presynaptically, Dlg4,
                 Homer1, Shank2, Shank3, Nlgn1, Camk2a postsynaptically. A
                 hand-made list is arguable, which is exactly why the next one
                 exists.
      psd_pc1    the first principal component of the 188 postsynaptic-density
                 genes (GO:0014069, less the subunit and localisation sets).
                 Nobody chose these individually: the component is whatever
                 those genes have most in common, and it is the stronger
                 density predictor of the two (0.26 against 0.15).
      autofluo   not a gene at all -- the autofluorescence of the same ten
                 brains, per structure. The only predictor measured in the same
                 tissue as the thing being predicted.

Outputs, under data\\adult_v2\\beyond\\for_sami:
  E_regression.png/.eps     the fit and its diagnostic
  F_maps.png/.eps           observed, predicted and residual on the brain
  regression_table.csv      every structure: observed, predicted, residual

  D:\\sep_histology\\code\\tools\\venv_atlas\\Scripts\\python.exe mapping\\v2_beyond_regression.py
"""

import csv
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import rankdata, spearmanr

from v2_per_mouse import annotation_20
from v2_beyond_density import (ADULTS, OUT, SUBUNITS, MARKERS, nano_per_mouse,
                               gene_profiles, autofluorescence, keep_structure,
                               build_covariates, half_map, residual, r_squared,
                               cv_r2, flexible, structure_names, tidy)
from v2_beyond_figures import RED, DARK, LIGHT, BLUE, FIGS, save

# Coronal planes to draw, in 20 um slices through the CCF (450 of them). Chosen
# to show cortex with hippocampus beneath it and thalamus at the midline, which
# is where the residual is largest in both directions.
PLANES = (215, 265, 315)

# The brain sits on black, as the young-versus-adult detail figures do. That
# matters more here than it looks: `hot` ends in white, so on a white page the
# brightest structures disappear into the background. It also BEGINS in black,
# though, so the floor is dropped below the data before drawing -- rank 1 then
# lands around a tenth of the way up the colormap, a dark red that still reads
# against the background instead of vanishing into it at the other end.
FLOOR = 0.12


def paint(plane, value_by_name, names):
    """One coronal slice with each structure filled by its value, NaN elsewhere."""
    labels = annotation_20('ccf')[plane]
    out = np.full(labels.shape, np.nan, np.float32)
    for idx in np.unique(labels):
        if idx == 0:
            continue
        v = value_by_name.get(names.get(int(idx), ''))
        if v is not None:
            out[labels == idx] = v
    return out


def panel_e(observed, predicted, res, structures, fitted_r2, cv):
    fig, axes = plt.subplots(1, 3, figsize=(13.4, 4.3))

    ax = axes[0]
    lim = [min(observed.min(), predicted.min()) - 4,
           max(observed.max(), predicted.max()) + 4]
    ax.plot(lim, lim, color=LIGHT, ls='--', lw=1.0, zorder=1)
    ax.scatter(predicted, observed, s=16, facecolor=DARK, edgecolor='0.2',
               linewidth=0.3, zorder=2)
    worst = np.argsort(-np.abs(res))[:5]
    for i in worst:
        ax.annotate(structures[i][:24], (predicted[i], observed[i]), fontsize=6.5,
                    color=RED, xytext=(4, 2), textcoords='offset points')
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel('predicted from abundance and density (rank)', fontsize=8.5)
    ax.set_ylabel('observed surface GluA1 (rank)', fontsize=8.5)
    ax.set_title(f'the fit\nR2 {fitted_r2:.2f} fitted, {cv:.2f} predicted', fontsize=9.5)
    tidy(ax)

    ax = axes[1]
    ax.axhline(0, color=LIGHT, lw=1.0)
    ax.scatter(predicted, res, s=16, facecolor=DARK, edgecolor='0.2',
               linewidth=0.3, zorder=2)
    ax.set_xlabel('predicted (rank)', fontsize=8.5)
    ax.set_ylabel('residual (ranks)', fontsize=8.5)
    ax.set_title('the diagnostic\nno tilt and no fan, so the model is\n'
                 'incomplete rather than mis-specified', fontsize=9.5)
    tidy(ax)

    ax = axes[2]
    ax.hist(res, bins=26, color=LIGHT, edgecolor='0.3', linewidth=0.4)
    ax.axvline(0, color=DARK, lw=1.4)
    ax.set_xlabel('residual (ranks)', fontsize=8.5)
    ax.set_ylabel('structures', fontsize=8.5)
    ax.set_title(f'the leftover\nspread {res.std():.1f} ranks over '
                 f'{len(res)} structures', fontsize=9.5)
    tidy(ax)

    fig.suptitle('E.  The regression behind the claim: surface GluA1 predicted from '
                 'receptor abundance and synaptic density', fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    save(fig, 'E_regression')


def panel_f(observed, predicted, res, structures):
    names = structure_names()
    maps = [('observed\nsurface GluA1', dict(zip(structures, observed)), 'hot', None),
            ('predicted from abundance\nand synaptic density',
             dict(zip(structures, predicted)), 'hot', None),
            ('what is left over',
             dict(zip(structures, res)), 'RdBu_r', float(np.abs(res).max()))]

    fig, axes = plt.subplots(len(PLANES), 3, figsize=(10.2, 3.2 * len(PLANES)))
    axes = np.atleast_2d(axes)
    for r, plane in enumerate(PLANES):
        for c, (title, values, cmap, span) in enumerate(maps):
            img = paint(plane, values, names)
            ax = axes[r, c]
            if span is None:
                lo = 1 - FLOOR * (len(structures) - 1)
                im = ax.imshow(img, cmap=cmap, vmin=lo, vmax=len(structures),
                               interpolation='nearest')
            else:
                im = ax.imshow(img, cmap=cmap, vmin=-span, vmax=span,
                               interpolation='nearest')
            im.set_rasterized(True)          # EPS keeps the text vector, not this
            ax.set_facecolor('black')        # NaN is transparent, so this is the ground
            ax.set_xticks([]); ax.set_yticks([])
            for side in ax.spines.values():
                side.set_visible(False)
            if r == 0:
                ax.set_title(title, fontsize=9)
            if c == 0:
                ax.set_ylabel(f'{plane * 0.02:.1f} mm', fontsize=8)
            if r == len(PLANES) - 1:
                cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02,
                                  orientation='horizontal')
                cb.ax.tick_params(labelsize=6.5)
                cb.set_label('rank among structures' if span is None
                             else 'observed minus predicted (ranks)', fontsize=7)
                if span is None:
                    # the floor sits below rank 1 so that nothing draws as black;
                    # the ticks should still stop at the real range
                    cb.set_ticks([1, 25, 50, 75, 100, len(structures)])

    fig.suptitle('F.  The same three quantities on the brain.  Red in the third column is '
                 'more surface GluA1 than\nabundance and density predict, blue is less.  '
                 'Black is outside the brain, or a structure the\nanalysis excludes: '
                 'fibre tracts, ventricles and unassigned voxels.', fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    save(fig, 'F_maps')


def main():
    os.makedirs(FIGS, exist_ok=True)
    nano, division = nano_per_mouse()
    expr, role = gene_profiles()
    everywhere = set.intersection(*[set(nano[m]) for m in ADULTS])
    everywhere &= set.intersection(*[set(expr[g]) for g in SUBUNITS + MARKERS])
    kept = [s for s in sorted(everywhere) if keep_structure(s, division.get(s, ''))[0]]
    auto = autofluorescence(set(kept))
    structures = sorted(set(kept).intersection(*[set(auto[m]) for m in ADULTS]))

    covariates, controls, _ = build_covariates(nano, expr, role, auto, structures)
    model = flexible(list(covariates.values()))
    observed = half_map(nano, range(len(ADULTS)), structures)
    res = residual(observed, model)
    predicted = observed - res
    fitted, cv = r_squared(observed, model), cv_r2(observed, model)

    print(f'{len(structures)} structures; predictors: Gria1-4 ({len(SUBUNITS)} genes), '
          f'{len(MARKERS)} hand-picked markers, psd_pc1 of {len(controls)} genes, '
          f'autofluorescence')
    print(f'  R2 {fitted:.3f} fitted, {cv:.3f} cross-validated; '
          f'observed against predicted rho {spearmanr(observed, predicted).statistic:.3f}')
    print(f'  residual spread {res.std():.1f} ranks; '
          f'residual against predicted rho {spearmanr(predicted, res).statistic:+.3f} '
          f'(should be ~0 if the model is not mis-specified)')

    with open(os.path.join(FIGS, 'regression_table.csv'), 'w', newline='',
              encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(['structure', 'observed_rank', 'predicted_rank', 'residual'])
        for i in np.argsort(-res):
            w.writerow([structures[i], f'{observed[i]:.1f}', f'{predicted[i]:.1f}',
                        f'{res[i]:.2f}'])

    panel_e(observed, predicted, res, structures, fitted, cv)
    panel_f(observed, predicted, res, structures)


if __name__ == '__main__':
    main()
