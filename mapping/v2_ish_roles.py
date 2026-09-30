"""
Subunit genes against localisation genes: the contrast that does discriminate.

The question is whether the nano map tracks GluA1 *at the membrane* or simply
GluA1. Two hypotheses:

  A   nano is proportional to total receptor
  B   nano is proportional to surface receptor = total x surface fraction

The postsynaptic-versus-presynaptic contrast in v2_ish_words.py cannot separate
these, because A and B both predict a postsynaptic map. But a contrast *within*
the postsynaptic compartment can, and this is it:

  what sets ABUNDANCE      the receptor's own subunit genes, Gria1-4
  what sets LOCALISATION   TARPs, cornichons, CKAMP44, the PSD scaffolds that
                           anchor the receptor, and NSF

Under A the subunit genes should carry the map and the localisation genes should
add nothing beyond them. Under B the map is a product of both, so the
localisation genes should explain variance that the subunits do not. The
asymmetry between those two statements is the measurement, and it is a category
statement, not a single-gene one -- which is the point: Cacng8 outranking Gria1
is one gene, and one gene can outrank another for any number of reasons.

**The curated roles below exist because the `category` column cannot express
this.** That column puts the metabotropic glutamate receptors Grm1-5 in
"auxiliary" beside the TARPs, and fills "trafficking" mostly with presynaptic
vesicle machinery -- Snap25, Stx1a, Bsn, Vamp2, Syt2 -- which is trafficking of
a completely different kind. The roles here are assigned by protein function,
from what each gene's product does, and no rho value was consulted in drawing
them. They were, however, written after the per-gene ranking had been seen, so
they are not blind; the permutation below is what guards the obvious version of
that worry.

Three things are computed, in increasing order of how much they are worth.

 1. **Rho by role.** The descriptive picture: how each curated set correlates
    with the map. This is the version that can be read off a figure.
 2. **Commonality.** Rank-regress the map on the subunit composite, on the
    localisation composite, and on both, and split the explained variance into
    what only the subunits explain, what only localisation explains, and what
    they share. The two composites are strongly co-expressed, so the shared
    part is large by construction -- the unique parts are the informative ones.
 3. **A within-family permutation.** The 19 AMPAR-family genes are split into 4
    and 15 in every possible way (3,876 of them) and the statistic recomputed.
    This asks whether the split by function is special *among splits of the same
    genes*, which controls for the co-expression that makes the conventional
    gene-category null useless (Fulcher 2021). It is exact, not sampled.

**PRESYN is the specificity control, and it matters.** If the nano map simply
likes genes for membrane trafficking of any kind, presynaptic vesicle machinery
should score too. If what it likes is AMPA receptor localisation specifically,
that set should not.

**What this still cannot do.** Gria1 mRNA is not GluA1 protein: translation,
transport and turnover vary regionally, so a subunit set that explains less than
localisation does is consistent with hypothesis B and also with a regional
translation gradient. Only a real total-receptor channel separates those, and
the one we had is not one (see v2_sep_channel_check.py). This is the strongest
version of the argument the current data supports, and it is suggestive.

Outputs, under data\\adult_v2\\ish:
  gene_roles.csv        the curated assignment, so it can be argued with
  role_summary.csv      per role and reading: n, median rho, quartiles
  ish_roles.png         the three panels above, for zref

  D:\\sep_histology\\code\\tools\\venv_atlas\\Scripts\\python.exe mapping\\v2_ish_roles.py
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

from v2_paths import DATA
NANO = os.path.join(DATA, 'comparisons_v2', 'young_vs_adult', 'region_means_per_mouse.csv')
GENES = os.path.join(DATA, 'adult_v2', 'ish', 'gene_region_table.csv')
OUT = os.path.join(DATA, 'adult_v2', 'ish')

ADULT_GROUPS = ('naive', 'rws')
READING = 'zref'                  # the reading the project has settled on
ALL_READINGS = ('zref', 'cref', 'subref', 'ratio', 'sepratio')
MIN_ISH_VOXELS = 10

# --------------------------------------------------------------- curated roles
# Assigned from what the protein does, not from any correlation. The two sets
# that carry the argument are SUBUNIT and LOCALISATION; PRESYN is the
# specificity control; the rest are context and are not part of any test.
ROLES = {
    # the receptor itself -- these set how much receptor there is
    'Gria1': 'subunit', 'Gria2': 'subunit', 'Gria3': 'subunit', 'Gria4': 'subunit',

    # auxiliary subunits: TARPs, cornichons, CKAMP44. They do not make receptor,
    # they decide how much of it reaches and stays in the synaptic membrane
    'Cacng3': 'localisation', 'Cacng4': 'localisation', 'Cacng5': 'localisation',
    'Cacng7': 'localisation', 'Cacng8': 'localisation',
    'Cnih2': 'localisation', 'Cnih3': 'localisation', 'Shisa9': 'localisation',
    # the PSD scaffolds that anchor AMPARs, and NSF/synaptojanin, which hold
    # them at the surface and retrieve them
    'Dlg4': 'localisation', 'Dlg2': 'localisation', 'Grip1': 'localisation',
    'Grip2': 'localisation', 'Homer1': 'localisation', 'Nsf': 'localisation',
    'Synj1': 'localisation',

    # vesicle fusion machinery: membrane trafficking, but presynaptic. The
    # control that separates "AMPAR localisation" from "trafficking in general"
    'Snap25': 'presyn', 'Stx1a': 'presyn', 'Syn1': 'presyn', 'Cplx1': 'presyn',
    'Bsn': 'presyn', 'Vamp2': 'presyn', 'Syt1': 'presyn', 'Syt2': 'presyn',
    'Syt3': 'presyn',

    # trans-synaptic organisers
    'Nrxn1': 'adhesion', 'Lrfn2': 'adhesion', 'Ptprs': 'adhesion',
    'Cbln2': 'adhesion', 'Nptx1': 'adhesion',

    # other glutamate receptors: NMDA, kainate, metabotropic
    'Grin1': 'other_glut_r', 'Grin2a': 'other_glut_r', 'Grin2b': 'other_glut_r',
    'Grin2c': 'other_glut_r', 'Grin2d': 'other_glut_r', 'Grin3a': 'other_glut_r',
    'Grik1': 'other_glut_r', 'Grik2': 'other_glut_r', 'Grik3': 'other_glut_r',
    'Grik4': 'other_glut_r', 'Grik5': 'other_glut_r',
    'Grm1': 'other_glut_r', 'Grm2': 'other_glut_r', 'Grm3': 'other_glut_r',
    'Grm4': 'other_glut_r', 'Grm5': 'other_glut_r',

    # excitatory identity, immediate-early genes, and the original controls
    'Slc17a7': 'excit_marker', 'Slc17a6': 'excit_marker', 'Camk2a': 'excit_marker',
    'Nrgn': 'excit_marker',
    'Arc': 'ieg', 'Bdnf': 'ieg', 'Egr1': 'ieg', 'Fos': 'ieg', 'Nos1': 'ieg',
}
FAMILY = ('subunit', 'localisation')     # the 19 genes the permutation shuffles
ORDER = ['subunit', 'localisation', 'presyn', 'adhesion', 'other_glut_r',
         'excit_marker', 'ieg', 'control_inhib', 'control_glia', 'control_struct']
NICE = {'subunit': 'AMPAR subunits\nGria1-4', 'localisation': 'AMPAR localisation\nTARPs, PSD, NSF',
        'presyn': 'presynaptic\nvesicle machinery', 'adhesion': 'trans-synaptic\norganisers',
        'other_glut_r': 'other glutamate\nreceptors', 'excit_marker': 'excitatory\nidentity',
        'ieg': 'immediate-early', 'control_inhib': 'inhibitory\ncontrols',
        'control_glia': 'glial\ncontrols', 'control_struct': 'structural\ncontrols'}


def adult_profile(reading):
    per = defaultdict(list)
    with open(NANO, newline='', encoding='utf-8') as fh:
        for r in csv.DictReader(fh):
            if r['group'] in ADULT_GROUPS and r['reading'] == reading:
                per[r['structure']].append(float(r['log2_value']))
    return {s: float(np.mean(v)) for s, v in per.items()}


def gene_profiles():
    out, cat = defaultdict(dict), {}
    with open(GENES, newline='', encoding='utf-8') as fh:
        for r in csv.DictReader(fh):
            if int(r['n_voxels']) >= MIN_ISH_VOXELS:
                out[r['symbol']][r['structure']] = float(r['ish_mean'])
                cat[r['symbol']] = r['category']
    return out, cat


def role_of(gene, category):
    """Curated role where there is one, the original control label otherwise."""
    if gene in ROLES:
        return ROLES[gene]
    c = category.get(gene, '')
    return c if c.startswith('control') else 'other'


def composite(genes, expr, structures):
    """One profile from a set of genes: the mean of their rank profiles."""
    r = [rankdata([expr[g][s] for s in structures]) for g in genes]
    return np.mean(r, axis=0)


def r2(y, xs):
    """Fraction of y explained by the columns in xs, linearly, on ranks."""
    a = np.column_stack(list(xs) + [np.ones_like(y)])
    resid = y - a @ np.linalg.lstsq(a, y, rcond=None)[0]
    return float(1 - resid.var() / y.var())


def commonality(y, s, m):
    """Unique to subunits, unique to localisation, and what they share."""
    rs, rm, both = r2(y, [s]), r2(y, [m]), r2(y, [s, m])
    return dict(r2_subunit=rs, r2_localisation=rm, r2_both=both,
                unique_subunit=both - rm, unique_localisation=both - rs,
                shared=rs + rm - both)


def permutation(y, expr, structures, family_genes, n_subunit, observed):
    """Every way of calling 4 of the 19 family genes the subunit set.

    The point of splitting the SAME genes is that co-expression, which is what
    makes the usual gene-category null meaningless here, is held fixed: each
    surrogate set is as co-expressed as the real one.
    """
    stats = []
    for pick in itertools.combinations(sorted(family_genes), n_subunit):
        rest = [g for g in family_genes if g not in pick]
        c = commonality(y, composite(pick, expr, structures),
                        composite(rest, expr, structures))
        stats.append(c['unique_localisation'] - c['unique_subunit'])
    stats = np.array(stats)
    p = float((np.sum(stats >= observed) + 1) / (len(stats) + 1))
    return stats, p


def sensitivity(y, expr, structures, family, roles):
    """The same test again, with the barely-expressed genes of the family gone.

    The obvious objection to a negative result here is that the localisation set
    is carrying passengers: Cacng5, Cacng7 and Grip2 are close to absent from the
    forebrain, so they can only dilute it. Dropping genes because they scored
    badly would be exactly the post-hoc move the permutation exists to catch, so
    the filter is on EXPRESSION ONLY -- median energy across structures, against
    the median of the family, applied to both sets alike and never looking at a
    rho. If the conclusion survives this, it is not an artefact of dead weight.
    """
    level = {g: float(np.median([expr[g][t] for t in structures])) for g in family}
    cut = float(np.median(list(level.values())))
    kept = sorted(g for g in family if level[g] >= cut)
    sub = [g for g in kept if roles[g] == 'subunit']
    loc = [g for g in kept if roles[g] == 'localisation']
    print(f'\n4. sensitivity: family genes at or above the family median expression '
          f'({cut:.1f})')
    print(f'   kept {len(kept)} of {len(family)}; dropped '
          f'{", ".join(g for g in family if g not in kept)}')
    if len(sub) < 2 or len(loc) < 2:
        print('   too few genes left on one side -- not run')
        return
    print(f'   subunit ({len(sub)}) {", ".join(sub)}')
    print(f'   localisation ({len(loc)}) {", ".join(loc)}')
    c = commonality(y, composite(sub, expr, structures), composite(loc, expr, structures))
    observed = c['unique_localisation'] - c['unique_subunit']
    stats, p = permutation(y, expr, structures, kept, len(sub), observed)
    print(f'   unique subunit {c["unique_subunit"]:+.3f}, '
          f'unique localisation {c["unique_localisation"]:+.3f}, '
          f'shared {c["shared"]:+.3f}')
    print(f'   observed {observed:+.4f} against {len(stats):,} splits, p = {p:.4f}')


def main():
    nano = adult_profile(READING)
    expr, category = gene_profiles()
    roles = {g: role_of(g, category) for g in expr}

    # every test uses the same structures: those the map and all 19 family
    # genes have in common, so nothing moves between comparisons
    family = sorted(g for g in expr if roles[g] in FAMILY)
    structures = sorted(set(nano).intersection(*[set(expr[g]) for g in family]))
    print(f'{len(expr)} genes, {len(family)} in the AMPAR family, '
          f'{len(structures)} shared structures, reading {READING}')
    sub = sorted(g for g in family if roles[g] == 'subunit')
    loc = sorted(g for g in family if roles[g] == 'localisation')
    print(f'  subunit      ({len(sub):2d})  {", ".join(sub)}')
    print(f'  localisation ({len(loc):2d})  {", ".join(loc)}')

    # ---------------------------------------------------------- 1. rho by role
    rows = []
    for reading in ALL_READINGS:
        prof = adult_profile(reading)
        for gene in sorted(expr):
            common = sorted(set(prof) & set(expr[gene]))
            if len(common) < 50:
                continue
            rho, _ = spearmanr([prof[s] for s in common], [expr[gene][s] for s in common])
            rows.append(dict(reading=reading, role=roles[gene], symbol=gene,
                             n_structures=len(common), rho=float(rho)))
    path = os.path.join(OUT, 'role_summary.csv')
    with open(path, 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows({k: (f'{v:.4f}' if isinstance(v, float) else v) for k, v in r.items()}
                    for r in rows)
    print(f'{len(rows)} rows -> {path}')

    with open(os.path.join(OUT, 'gene_roles.csv'), 'w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh); w.writerow(['symbol', 'role', 'original_category'])
        for g in sorted(expr):
            w.writerow([g, roles[g], category.get(g, '')])

    here = [r for r in rows if r['reading'] == READING]
    by_role = defaultdict(list)
    for r in here:
        by_role[r['role']].append(r['rho'])
    print(f'\n1. rho with the {READING} map, by curated role')
    for role in ORDER:
        v = by_role.get(role, [])
        if v:
            print(f'  {role:15s} n={len(v):3d}   median {np.median(v):+.3f}   '
                  f'[{np.percentile(v, 25):+.3f} {np.percentile(v, 75):+.3f}]')

    # -------------------------------------------------------- 2. commonality
    y = rankdata([nano[s] for s in structures])
    s_comp = composite(sub, expr, structures)
    m_comp = composite(loc, expr, structures)
    c = commonality(y, s_comp, m_comp)
    print(f'\n2. variance of the {READING} map explained (ranks, {len(structures)} structures)')
    print(f'  subunit composite alone        R2 = {c["r2_subunit"]:.3f}')
    print(f'  localisation composite alone   R2 = {c["r2_localisation"]:.3f}')
    print(f'  both                           R2 = {c["r2_both"]:.3f}')
    print(f'  unique to subunits             {c["unique_subunit"]:+.3f}')
    print(f'  unique to localisation         {c["unique_localisation"]:+.3f}')
    print(f'  shared                         {c["shared"]:+.3f}')
    print(f'  composites correlate with each other at rho '
          f'{spearmanr(s_comp, m_comp).statistic:+.3f}')

    # -------------------------------------------------------- 3. permutation
    observed = c['unique_localisation'] - c['unique_subunit']
    stats, p = permutation(y, expr, structures, family, len(sub), observed)
    print(f'\n3. every 4-of-{len(family)} split of the same family: {len(stats):,} of them')
    print(f'  observed difference in unique variance {observed:+.4f}')
    print(f'  surrogate splits median {np.median(stats):+.4f}, '
          f'p95 {np.percentile(stats, 95):+.4f}')
    print(f'  p = {p:.4f}  (fraction of splits reaching the functional one)')

    # the specificity control
    pre = by_role.get('presyn', [])
    print(f'\n  specificity control -- presynaptic vesicle machinery: '
          f'n={len(pre)}, median rho {np.median(pre):+.3f} '
          f'against localisation {np.median(by_role["localisation"]):+.3f}')

    sensitivity(y, expr, structures, family, roles)
    figure(by_role, c, stats, observed, p, s_comp, m_comp, y)


def figure(by_role, c, stats, observed, p, s_comp, m_comp, y):
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.2),
                             gridspec_kw=dict(width_ratios=[1.7, 0.8, 1.0]))
    rng = np.random.default_rng(0)

    ax = axes[0]
    roles = [r for r in ORDER if by_role.get(r)]
    for i, role in enumerate(roles):
        v = by_role[role]
        colour = {'subunit': '#1f3b73', 'localisation': '#c0392b'}.get(role, '0.65')
        ax.scatter(np.full(len(v), i) + rng.uniform(-.14, .14, len(v)), v, s=17,
                   facecolor=colour, edgecolor='0.25', linewidth=0.4, zorder=2)
        ax.plot([i - .3, i + .3], [np.median(v)] * 2, color='0.15', lw=1.8, zorder=3)
    ax.axhline(0, color='0.85', lw=0.7, zorder=0)
    ax.set_xticks(range(len(roles)))
    ax.set_xticklabels([f'{NICE.get(r, r).replace(chr(10), " ")} ({len(by_role[r])})'
                        for r in roles], fontsize=7, rotation=32, ha='right')
    ax.set_ylabel(f'Spearman with the adult {READING} map', fontsize=8)
    ax.set_title('1. the map by curated role\n'
                 'blue = what sets abundance, red = what sets localisation', fontsize=9)

    ax = axes[1]
    parts = [c['unique_subunit'], c['shared'], c['unique_localisation']]
    ax.bar(range(3), parts, color=['#1f3b73', '0.75', '#c0392b'],
           edgecolor='0.25', linewidth=0.5)
    ax.set_xticks(range(3))
    ax.set_xticklabels(['unique to\nsubunits', 'shared', 'unique to\nlocalisation'],
                       fontsize=7.5)
    ax.set_ylabel('variance of the map explained (R2, ranks)', fontsize=8)
    for i, v in enumerate(parts):
        ax.text(i, v, f'{v:.3f}', ha='center', va='bottom', fontsize=7.5)
    ax.set_title(f'2. commonality\ncomposites correlate at rho '
                 f'{spearmanr(s_comp, m_comp).statistic:+.2f}', fontsize=9)

    ax = axes[2]
    ax.hist(stats, bins=60, color='0.72', edgecolor='0.35', linewidth=0.3)
    ax.axvline(observed, color='#c0392b', lw=2)
    ax.annotate('the split\nby function', (observed, ax.get_ylim()[1] * 0.9),
                color='#c0392b', fontsize=7.5, ha='right', va='top',
                xytext=(-6, 0), textcoords='offset points')
    ax.set_xlabel('unique(localisation) - unique(subunit)', fontsize=8)
    ax.set_ylabel(f'splits of the same {len(stats):,}', fontsize=8)
    ax.set_title(f'3. every 4-of-19 split of the family\np = {p:.4f}', fontsize=9)

    for ax in axes:
        ax.tick_params(labelsize=7)
        for side in ('top', 'right'):
            ax.spines[side].set_visible(False)
    fig.suptitle('Does the map need the localisation genes, or do the subunits account for it?'
                 '\nThe direction is as predicted -- the subunits add almost nothing beyond '
                 'localisation -- but splitting the same 19 genes any other way does as well '
                 '(panel 3), so this panel is no evidence either way.', fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.89))
    path = os.path.join(OUT, 'ish_roles.png')
    fig.savefig(path, dpi=200); plt.close(fig)
    print(f'\n{path}')


if __name__ == '__main__':
    main()
