"""
The gene panel, defined by ontology rather than by us.

The subunit-against-localisation test in v2_ish_roles.py came out with the
predicted direction and no evidence for it, because the panel is too small: four
subunit genes, fifteen localisation genes, and a permutation null wide enough to
swallow the effect four times over. This builds the panel that could carry it.

**Nothing here is hand-picked.** Membership comes from Gene Ontology terms,
named below with the reason each is in the list, queried through mygene.info and
cached. Every gene in the output carries the term ids that put it there, so the
panel can be argued with term by term rather than gene by gene -- which is the
whole point, since a panel we curate ourselves is a panel that can agree with us.

  subunit       what sets how much receptor there is
                GO:0004971 AMPA glutamate receptor activity, intersected with
                GO:0032281 AMPA glutamate receptor complex, then MINUS two genes
                by hand -- see OVERRIDE below, the only hand-made call in this
                file and the only one in the panel.
  localisation  what decides where the receptor goes and whether it stays
                GO:0099072 regulation of postsynaptic membrane neurotransmitter
                           receptor levels
                GO:0099645 neurotransmitter receptor localization to
                           postsynaptic specialization membrane
                GO:0097113 AMPA glutamate receptor clustering
                GO:0098970 postsynaptic neurotransmitter receptor diffusion
                           trapping
                GO:0032281 AMPA glutamate receptor complex, minus the subunits
                           -- what is left is the auxiliary subunits, the TARPs
                           and cornichons
  control_psd   the matched control pool: GO:0014069 postsynaptic density,
                minus everything above. These are postsynaptic genes with no
                annotated role in getting receptors to the membrane, which is
                exactly the comparison the test needs -- it is what separates
                "the map likes postsynaptic genes" from "the map likes receptor
                localisation genes".

A gene claimed by more than one set goes to the most specific: subunit, then
localisation, then control.

**Both planes of section are kept, and that is a change.** The old panel was
coronal only. Most Allen genes have only a sagittal experiment, so coronal-only
would throw away most of the panel; sagittal grids sit in the same reference
space and cover one hemisphere, which is enough for a structure mean that pools
both. The reason to want them is stronger than that, though: a gene with a
coronal AND a sagittal experiment has been measured twice, independently, so
those genes give a direct estimate of how reliable one Allen map is -- the
question standing open since the first version of this analysis. One row per
experiment here, not per gene, so that stays visible downstream.

Outputs, under data\\adult_v2\\panel:
  panel_v2.csv        one row per experiment: symbol, role, go_terms, id, plane
  panel_genes.csv     one row per gene, with every term that claimed it
  cache\\*.json        every API response, so a re-run asks nothing twice

  D:\\sep_histology\\code\\tools\\venv_atlas\\Scripts\\python.exe v2_panel_build.py
"""

import csv
import json
import os
import time
import urllib.parse
import urllib.request
from collections import defaultdict

from v2_paths import DATA
OUT = os.path.join(DATA, 'adult_v2', 'panel')
CACHE = os.path.join(OUT, 'cache')
OLD_PANEL = os.path.join(DATA, 'gene_targets.csv')

MYGENE = 'https://mygene.info/v3/query'
ALLEN = 'http://api.brain-map.org/api/v2/data/query.json'

# term -> (role it contributes to, what it means). Order matters only for
# readability; precedence is set by ROLE_ORDER below.
TERMS = {
    'GO:0004971': ('subunit_candidate', 'AMPA glutamate receptor activity'),
    'GO:0032281': ('ampar_complex', 'AMPA glutamate receptor complex'),
    'GO:0099072': ('localisation', 'regulation of postsynaptic membrane neurotransmitter receptor levels'),
    'GO:0099645': ('localisation', 'neurotransmitter receptor localization to postsynaptic specialization membrane'),
    'GO:0097113': ('localisation', 'AMPA glutamate receptor clustering'),
    'GO:0098970': ('localisation', 'postsynaptic neurotransmitter receptor diffusion trapping'),
    'GO:0014069': ('control_psd', 'postsynaptic density'),
}
# THE ONE HAND-MADE DECISION IN THIS FILE, written down so it can be argued
# with and measured. GO annotates the delta-family receptors Grid1 and Grid2 to
# both AMPA receptor terms, so the intersection above picks them up. They are a
# separate family of ionotropic receptors: they do not form AMPA receptors and
# their transcript level does not set how much AMPA receptor a region has, which
# is the only property the subunit set is supposed to carry. They are kept in the
# panel under their own role rather than dropped, so the test can be run with
# them counted as subunits and the difference reported.
OVERRIDE = {'Grid1': 'delta_receptor', 'Grid2': 'delta_receptor'}

ROLE_ORDER = ('subunit', 'delta_receptor', 'localisation', 'control_psd')
MAX_HITS = 1000


def cached(name, fetch):
    """Every API answer lands on disk before it is used."""
    os.makedirs(CACHE, exist_ok=True)
    path = os.path.join(CACHE, name + '.json')
    if os.path.exists(path):
        with open(path, encoding='utf-8') as fh:
            return json.load(fh)
    value = fetch()
    with open(path, 'w', encoding='utf-8') as fh:
        json.dump(value, fh, indent=1)
    time.sleep(0.2)
    return value


def genes_with_term(term):
    """Every mouse gene annotated to one GO term."""
    def fetch():
        query = urllib.parse.urlencode({'q': 'go:' + term.split(':')[-1], 'species': 'mouse',
                                        'fields': 'symbol', 'size': MAX_HITS})
        with urllib.request.urlopen(f'{MYGENE}?{query}', timeout=90) as fh:
            d = json.load(fh)
        if d.get('total', 0) > MAX_HITS:
            raise SystemExit(f'{term}: {d["total"]} genes, above the {MAX_HITS} cap')
        return sorted({h['symbol'] for h in d.get('hits', []) if h.get('symbol')})
    return cached('term_' + term.replace(':', '_'), fetch)


def allen_experiments(symbol):
    """Every usable Allen ISH experiment for one gene, with its plane."""
    def fetch():
        crit = ("model::SectionDataSet,"
                "rma::criteria,[failed$eq'false'],products[abbreviation$eq'Mouse'],"
                f"genes[acronym$eq'{symbol}'],"
                "rma::include,plane_of_section,genes")
        query = urllib.parse.urlencode({'criteria': crit, 'num_rows': 50})
        with urllib.request.urlopen(f'{ALLEN}?{query}', timeout=90) as fh:
            d = json.load(fh)
        return [{'id': r['id'], 'plane': r['plane_of_section']['name']}
                for r in d.get('msg', []) if r.get('plane_of_section')]
    return cached('allen_' + symbol, fetch)


def assign_roles():
    """{gene: (role, [terms that claimed it])}, from the ontology alone."""
    members = {term: set(genes_with_term(term)) for term in TERMS}
    for term, genes in members.items():
        print(f'  {term}  {len(genes):4d} genes  {TERMS[term][1]}')

    subunit = members['GO:0004971'] & members['GO:0032281']
    localisation = set().union(*[members[t] for t in TERMS
                                 if TERMS[t][0] == 'localisation'])
    localisation |= members['GO:0032281'] - subunit      # the auxiliary subunits
    localisation -= subunit
    control = members['GO:0014069'] - subunit - localisation

    roles, why = {}, defaultdict(list)
    for role, genes in (('subunit', subunit), ('localisation', localisation),
                        ('control_psd', control)):
        for g in genes:
            roles[g] = role
    for g, role in OVERRIDE.items():
        if g in roles:
            print(f'  override: {g} {roles[g]} -> {role}')
            roles[g] = role
    for term, genes in members.items():
        for g in genes:
            if g in roles:
                why[g].append(term)
    return roles, why


def main():
    os.makedirs(OUT, exist_ok=True)
    print('gene sets, straight from the ontology:')
    roles, why = assign_roles()
    counts = defaultdict(int)
    for r in roles.values():
        counts[r] += 1
    print('\nassigned:', ', '.join(f'{r} {counts[r]}' for r in ROLE_ORDER))
    print('  subunit set is', ', '.join(sorted(g for g in roles if roles[g] == 'subunit')))

    # keep whatever the old panel already paid for, so nothing is re-downloaded
    old = {}
    with open(OLD_PANEL, newline='', encoding='utf-8') as fh:
        for r in csv.DictReader(fh):
            old[r['symbol']] = r['category']

    rows, per_gene, no_exp = [], [], []
    for i, gene in enumerate(sorted(roles), 1):
        exps = allen_experiments(gene)
        if not exps:
            no_exp.append(gene)
            continue
        planes = sorted({e['plane'] for e in exps})
        per_gene.append(dict(symbol=gene, role=roles[gene],
                             go_terms=' '.join(sorted(why[gene])),
                             n_experiments=len(exps), planes=' '.join(planes),
                             in_old_panel=old.get(gene, '')))
        for e in exps:
            rows.append(dict(symbol=gene, role=roles[gene], experiment_id=e['id'],
                             plane=e['plane'], category=roles[gene],
                             go_terms=' '.join(sorted(why[gene])),
                             in_old_panel=old.get(gene, '')))
        if i % 25 == 0:
            print(f'  {i}/{len(roles)} genes resolved', flush=True)

    path = os.path.join(OUT, 'panel_v2.csv')
    with open(path, 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    with open(os.path.join(OUT, 'panel_genes.csv'), 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=list(per_gene[0].keys()))
        w.writeheader(); w.writerows(per_gene)

    both = [g for g in per_gene if len(g['planes'].split()) > 1]
    print(f'\n{len(per_gene)} genes with at least one experiment, '
          f'{len(rows)} experiments -> {path}')
    for role in ROLE_ORDER:
        sel = [g for g in per_gene if g['role'] == role]
        exp = sum(g['n_experiments'] for g in sel)
        print(f'  {role:14s} {len(sel):4d} genes, {exp:4d} experiments')
    print(f'  {len(both)} genes measured in both planes -- these are what the '
          f'reliability estimate rests on')
    print(f'  {len(no_exp)} genes have no usable experiment'
          + (f': {", ".join(no_exp[:12])}...' if len(no_exp) > 12 else
             (f': {", ".join(no_exp)}' if no_exp else '')))
    already = sum(1 for r in rows if os.path.exists(
        os.path.join(DATA, 'atlas_ish', f'{r["experiment_id"]}_energy.mhd')))
    print(f'  {already} of {len(rows)} grids are already on disk')


if __name__ == '__main__':
    main()
