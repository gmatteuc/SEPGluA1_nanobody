"""The kinds of genes of analysis 3, fixed before any correlation with the map is read.

The April result grouped genes by categories written while the genes were being
chosen (gene_targets.csv), and split one of them by hand after looking; its
p = 0.032 rested on that split. Here the groups come from annotation nobody here
wrote, or from a cited marker list, and this file was committed before any gene
was correlated with the map on the new inputs (git log shows the order):

    subunits             GO:0004971 AMPA glutamate receptor activity, intersected
                         with GO:0032281 AMPA glutamate receptor complex, minus
                         Grid1 and Grid2 (ish.panel_build's rule, role subunit)
    localisation         GO:0099072, GO:0099645, GO:0097113, GO:0098970, and
                         GO:0032281 minus the subunits (panel_build's rule, role
                         localisation): the genes that bring AMPA receptors to the
                         postsynaptic membrane and hold them there
    other postsynaptic   annotated (cellular component) to GO:0098794 postsynapse
                         or any term below it, by is_a or part_of, and not to the
                         presynapse; not a subunit or localisation gene
    presynaptic          annotated to GO:0098793 presynapse or below, not to the
                         postsynapse; not a subunit or localisation gene, so the
                         four GO sets share no gene
    GABAergic markers    GABA synthesis and vesicular transport (Gad1, Gad2,
                         Slc32a1), and the GABAergic subclass markers of Tasic et
                         al. 2018 that have a grid here (of Lamp5, Sncg, Serpinf1,
                         Vip, Sst, Pvalb and Meis2: Sst and Pvalb)
    glia                 P9's control_glia, unchanged (Aldh1l1, Aqp4, Gfap, Olig2,
                         Sox10, Mbp, Cx3cr1)

GO is read with its ancestry, from a dated go-basic.obo (ish.gene_table), so a gene
annotated to "postsynaptic density" counts as postsynaptic. Annotations with a NOT
qualifier are left out. SynGO's curated annotations enter GO, so the GO sets carry
those that reached the records of mygene.info; SynGO's own release is not read (it
ships as spreadsheets, and venv_atlas has no reader for them). GO has no cell-type
terms, so the two control sets come from a cited list and from P9's own glia list,
neither chosen by correlation. A gene may sit in one GO set and one marker set
(Slc32a1 is presynaptic by GO and a GABAergic marker); a contrast between two sets
leaves out the genes in both, and says so.

The tests of the sets against the map come with analysis 3 (run_ish_gene_sets).
Membership is computed by run_ish_gene_table.py, with the gene table.
"""

# the sets in the order the figures draw them
SET_ORDER = (
    "subunits",
    "localisation",
    "other postsynaptic",
    "presynaptic",
    "GABAergic markers",
    "glia",
)

# the roles of ish.panel_build that define the first two sets
PANEL_ROLE_SETS = {"subunits": "subunit", "localisation": "localisation"}

# the two cellular components of the next two sets
POSTSYNAPSE = "GO:0098794"
PRESYNAPSE = "GO:0098793"

# GABA synthesis and vesicular transport, and the GABAergic subclass markers named by
# Tasic et al. 2018 (Nature 563:72); only those with a grid enter the set
GABA_MACHINERY = ("Gad1", "Gad2", "Slc32a1")
TASIC_GABA_MARKERS = ("Lamp5", "Sncg", "Serpinf1", "Vip", "Sst", "Pvalb", "Meis2")

# P9's control_glia, unchanged
GLIA = ("Aldh1l1", "Aqp4", "Gfap", "Olig2", "Sox10", "Mbp", "Cx3cr1")

# where each set comes from, as the gene table and the figures write it
GENE_SETS = {
    "subunits": "GO:0004971 and GO:0032281, minus Grid1 and Grid2 (panel role subunit)",
    "localisation": (
        "GO:0099072, GO:0099645, GO:0097113, GO:0098970, GO:0032281 minus the "
        "subunits (panel role localisation)"
    ),
    "other postsynaptic": (
        "GO:0098794 postsynapse or below, not presynapse, not in the two sets above"
    ),
    "presynaptic": (
        "GO:0098793 presynapse or below, not postsynapse, not in the first two sets"
    ),
    "GABAergic markers": "Gad1, Gad2, Slc32a1, and Tasic et al. 2018 subclass markers",
    "glia": "P9's control_glia: Aldh1l1, Aqp4, Gfap, Olig2, Sox10, Mbp, Cx3cr1",
}


def gene_sets(
    genes: list[str], role: dict[str, str], components: dict[str, set[str]]
) -> dict[str, list[str]]:
    """The members of each set among `genes`, {set: sorted symbols}.

    `role` is each gene's role in the ontology panel ('' for a gene outside it),
    `components` each gene's cellular-component terms with their ancestors (empty
    for a gene without a GO record).
    """
    out = {}
    for name, panel_role in PANEL_ROLE_SETS.items():
        out[name] = sorted(g for g in genes if role.get(g, "") == panel_role)
    by_role = set(out["subunits"]) | set(out["localisation"])
    post = {g for g in genes if POSTSYNAPSE in components.get(g, set())}
    pre = {g for g in genes if PRESYNAPSE in components.get(g, set())}
    out["other postsynaptic"] = sorted(post - pre - by_role)
    out["presynaptic"] = sorted(pre - post - by_role)
    markers = GABA_MACHINERY + TASIC_GABA_MARKERS
    out["GABAergic markers"] = sorted(g for g in genes if g in markers)
    out["glia"] = sorted(g for g in genes if g in GLIA)
    return {name: out[name] for name in SET_ORDER}
