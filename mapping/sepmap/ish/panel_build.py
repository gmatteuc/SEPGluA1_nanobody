"""The ontology-defined ISH gene panel: genes from GO terms, experiments from Allen.

The subunit-against-localisation test of ish.roles had too small a panel: four
subunit genes, fifteen localisation genes, and a permutation null four times wider
than the effect. This builds a panel that could carry it, with no gene picked by
hand: membership comes from Gene Ontology terms (TERMS, each with the reason it is
in the list), queried through mygene.info and cached, and every gene carries the
terms that put it there, so the panel can be argued with term by term. A panel
curated by us could agree with us.

    subunit         GO:0004971 AMPA glutamate receptor activity, intersected with
                    GO:0032281 AMPA glutamate receptor complex, minus Grid1 and
                    Grid2 (OVERRIDE, the one hand-made call of the panel)
    localisation    GO:0099072, GO:0099645, GO:0097113 and GO:0098970 (receptor
                    levels, localisation, clustering and diffusion trapping at the
                    postsynaptic membrane), and GO:0032281 minus the subunits,
                    which leaves the auxiliary subunits, TARPs and cornichons
    control_psd     GO:0014069 postsynaptic density minus everything above:
                    postsynaptic genes with no annotated role in getting receptors
                    to the membrane, which separates "the map likes postsynaptic
                    genes" from "the map likes receptor localisation genes"

A gene claimed by more than one set goes to the most specific: subunit, then
localisation, then control. Both planes of section are kept: most Allen genes have
only a sagittal experiment, and a sagittal grid sits in the same reference space
and covers one hemisphere, enough for a structure mean that pools both. A gene with
a coronal and a sagittal experiment is measured twice, independently, which gives
the reliability of one Allen map (ish.reliability), so the panel has one row per
experiment, not per gene.

Run by run_panel_build.py.
"""

import csv
import json
import os
import time
import urllib.parse
import urllib.request
from collections import defaultdict

from sepmap.config import DATA

OUT = os.path.join(DATA, "adult_v2", "panel")
CACHE = os.path.join(OUT, "cache")

# the 100-gene panel, whose category of each gene is recorded beside its new role
OLD_PANEL = os.path.join(DATA, "gene_targets.csv")

MYGENE = "https://mygene.info/v3/query"
ALLEN = "http://api.brain-map.org/api/v2/data/query.json"

# GO term: (the set it contributes to, what it means); assign_roles decides which
# set wins for a gene claimed by several
TERMS = {
    "GO:0004971": ("subunit_candidate", "AMPA glutamate receptor activity"),
    "GO:0032281": ("ampar_complex", "AMPA glutamate receptor complex"),
    "GO:0099072": (
        "localisation",
        "regulation of postsynaptic membrane neurotransmitter receptor levels",
    ),
    "GO:0099645": (
        "localisation",
        "neurotransmitter receptor localization to postsynaptic specialization membrane",
    ),
    "GO:0097113": ("localisation", "AMPA glutamate receptor clustering"),
    "GO:0098970": (
        "localisation",
        "postsynaptic neurotransmitter receptor diffusion trapping",
    ),
    "GO:0014069": ("control_psd", "postsynaptic density"),
}

# the one hand-made call of the panel, written down so it can be argued with: GO
# annotates the delta-family receptors Grid1 and Grid2 to both AMPA receptor terms,
# so the intersection picks them up, but they are a separate family of ionotropic
# receptors that do not form AMPA receptors, and their transcript level does not set
# how much AMPA receptor a region has, the one property the subunit set is to carry;
# they keep a role of their own rather than being dropped, so the test can also be
# run with them counted as subunits
OVERRIDE = {"Grid1": "delta_receptor", "Grid2": "delta_receptor"}

# the roles in the order they are printed
ROLE_ORDER = ("subunit", "delta_receptor", "localisation", "control_psd")

# genes asked of mygene.info per term; a term with more stops the run
MAX_HITS = 1000


def cached(name, fetch):
    """The answer of `fetch()`, cached as <name>.json, so a re-run asks nothing twice."""
    os.makedirs(CACHE, exist_ok=True)
    path = os.path.join(CACHE, name + ".json")
    if os.path.exists(path):
        with open(path, encoding="utf-8") as fh:
            return json.load(fh)
    value = fetch()
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(value, fh, indent=1)

    # a courtesy to the public APIs
    time.sleep(0.2)
    return value


def genes_with_term(term):
    """Every mouse gene annotated to one GO term, as a sorted list of symbols."""

    def fetch():
        """Ask mygene.info for the term's genes."""
        query = urllib.parse.urlencode(
            {
                "q": "go:" + term.split(":")[-1],
                "species": "mouse",
                "fields": "symbol",
                "size": MAX_HITS,
            }
        )
        with urllib.request.urlopen(f"{MYGENE}?{query}", timeout=90) as fh:
            d = json.load(fh)
        if d.get("total", 0) > MAX_HITS:
            raise SystemExit(f"{term}: {d['total']} genes, above the {MAX_HITS} cap")
        return sorted({h["symbol"] for h in d.get("hits", []) if h.get("symbol")})

    return cached("term_" + term.replace(":", "_"), fetch)


def allen_experiments(symbol):
    """Every usable Allen ISH experiment of one gene: [{"id": ..., "plane": ...}].

    Usable means not failed, from the mouse product, with a plane of section.
    """

    def fetch():
        """Ask the Allen API for the gene's experiments."""
        crit = (
            "model::SectionDataSet,"
            "rma::criteria,[failed$eq'false'],products[abbreviation$eq'Mouse'],"
            f"genes[acronym$eq'{symbol}'],"
            "rma::include,plane_of_section,genes"
        )
        query = urllib.parse.urlencode({"criteria": crit, "num_rows": 50})
        with urllib.request.urlopen(f"{ALLEN}?{query}", timeout=90) as fh:
            d = json.load(fh)
        return [
            {"id": r["id"], "plane": r["plane_of_section"]["name"]}
            for r in d.get("msg", [])
            if r.get("plane_of_section")
        ]

    return cached("allen_" + symbol, fetch)


def assign_roles():
    """{gene: role} and {gene: [terms that claimed it]}, from the ontology alone."""
    # the genes of each term
    members = {term: set(genes_with_term(term)) for term in TERMS}
    for term, genes in members.items():
        print(f"  {term}  {len(genes):4d} genes  {TERMS[term][1]}")

    # the three sets, the most specific first; GO:0032281 minus the subunits is the
    # auxiliary subunits
    subunit = members["GO:0004971"] & members["GO:0032281"]
    localisation = set().union(
        *[members[t] for t in TERMS if TERMS[t][0] == "localisation"]
    )
    localisation |= members["GO:0032281"] - subunit
    localisation -= subunit
    control = members["GO:0014069"] - subunit - localisation

    # each gene's role, the hand-made call, and the terms behind each gene
    roles, why = {}, defaultdict(list)
    for role, genes in (
        ("subunit", subunit),
        ("localisation", localisation),
        ("control_psd", control),
    ):
        for g in genes:
            roles[g] = role
    for g, role in OVERRIDE.items():
        if g in roles:
            print(f"  override: {g} {roles[g]} -> {role}")
            roles[g] = role
    for term, genes in members.items():
        for g in genes:
            if g in roles:
                why[g].append(term)
    return roles, why


def main():
    """Build the panel from the ontology and the Allen API; write both tables."""
    # the gene sets
    os.makedirs(OUT, exist_ok=True)
    print("gene sets, straight from the ontology:")
    roles, why = assign_roles()
    counts = defaultdict(int)
    for r in roles.values():
        counts[r] += 1
    print("\nassigned:", ", ".join(f"{r} {counts[r]}" for r in ROLE_ORDER))
    print(
        "  subunit set is", ", ".join(sorted(g for g in roles if roles[g] == "subunit"))
    )

    # the old panel's category of each gene, recorded beside the new role
    old = {}
    with open(OLD_PANEL, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            old[r["symbol"]] = r["category"]

    # the Allen experiments of each gene
    rows, per_gene, no_exp = [], [], []
    for i, gene in enumerate(sorted(roles), 1):
        exps = allen_experiments(gene)
        if not exps:
            no_exp.append(gene)
            continue
        planes = sorted({e["plane"] for e in exps})
        per_gene.append(
            dict(
                symbol=gene,
                role=roles[gene],
                go_terms=" ".join(sorted(why[gene])),
                n_experiments=len(exps),
                planes=" ".join(planes),
                in_old_panel=old.get(gene, ""),
            )
        )
        for e in exps:
            rows.append(
                dict(
                    symbol=gene,
                    role=roles[gene],
                    experiment_id=e["id"],
                    plane=e["plane"],
                    category=roles[gene],
                    go_terms=" ".join(sorted(why[gene])),
                    in_old_panel=old.get(gene, ""),
                )
            )
        if i % 25 == 0:
            print(f"  {i}/{len(roles)} genes resolved", flush=True)

    # one row per experiment, one row per gene
    path = os.path.join(OUT, "panel_v2.csv")
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    with open(
        os.path.join(OUT, "panel_genes.csv"), "w", newline="", encoding="utf-8"
    ) as fh:
        w = csv.DictWriter(fh, fieldnames=list(per_gene[0].keys()))
        w.writeheader()
        w.writerows(per_gene)

    # what the panel holds, and how much of it is already downloaded
    both = [g for g in per_gene if len(g["planes"].split()) > 1]
    print(
        f"\n{len(per_gene)} genes with at least one experiment, "
        f"{len(rows)} experiments -> {path}"
    )
    for role in ROLE_ORDER:
        sel = [g for g in per_gene if g["role"] == role]
        exp = sum(g["n_experiments"] for g in sel)
        print(f"  {role:14s} {len(sel):4d} genes, {exp:4d} experiments")
    print(
        f"  {len(both)} genes measured in both planes -- these are what the "
        f"reliability estimate rests on"
    )
    print(
        f"  {len(no_exp)} genes have no usable experiment"
        + (
            f": {', '.join(no_exp[:12])}..."
            if len(no_exp) > 12
            else (f": {', '.join(no_exp)}" if no_exp else "")
        )
    )
    already = sum(
        1
        for r in rows
        if os.path.exists(
            os.path.join(DATA, "atlas_ish", f"{r['experiment_id']}_energy.mhd")
        )
    )
    print(f"  {already} of {len(rows)} grids are already on disk")
