"""The guided figures of the ISH line: their numbers, questions and file names.

The figures are numbered in the order of the argument: the question (00), the inputs
(01, 02), part 1, how much of the map Gria1 and synapse density leave (03, 04), part
2, the genes against the map and against the leftover (05 to 11), the controls (12
to 14), the limit (15) and April's headline (16). A main figure's detailed version
carries its number with an s. Every figure's title, its file name and the figures'
references to each other read the tables here, and so do the index of the figures
(figures/README.md) and the overview figure that ish.overview writes: the parts of
the walk, the detailed versions of each main figure and the run script that draws
each. The docstrings and the run scripts' headers that name the files follow them by
hand.

Imported by the modules that draw a guided figure and by ish.overview.
"""

from pathlib import Path

from sepmap.structures import FIGURES

# the number of each guided figure, by its key
FIGURE_NUMBERS = {
    "overview": "00",
    "structures": "01",
    "genes": "02",
    "genes_detail": "02s",
    "beyond": "03",
    "beyond_budget": "03s1",
    "beyond_controls": "03s2",
    "beyond_where": "04",
    "one_comparison": "05",
    "one_comparison_detail": "05s",
    "spatial_null": "06",
    "spatial_null_detail": "06s",
    "top_genes": "07",
    "gene_ranking": "07s",
    "cacng8_gria1": "08",
    "between_within": "09",
    "between_within_detail": "09s",
    "gene_kinds": "10",
    "gene_sets": "10s1",
    "localisation": "10s2",
    "leftover": "11",
    "leftover_genes": "11s1",
    "ampa_family": "11s2",
    "autofluorescence": "12",
    "autofluorescence_detail": "12s",
    "robustness": "13",
    "robustness_detail": "13s",
    "synaptome": "14",
    "synaptome_detail": "14s1",
    "density_markers": "14s2",
    "green_channel": "15",
    "green_channel_detail": "15s",
    "april_headline": "16",
    "april_headline_detail": "16s",
}

# the question each guided figure answers: its title, and its heading in the index of
# the figures (ish.overview)
QUESTIONS = {
    "overview": "The ISH line on one page",
    "structures": "Which structures enter every comparison?",
    "genes": "Which genes, and how good is one Allen map?",
    "genes_detail": "The genes in detail: experiments, section QC, reliability, what "
    "was left out",
    "beyond": "Do Gria1 expression and synapse density explain the map?",
    "beyond_budget": "Part 1 in detail: each term, the partition, the calibration, "
    "the check rows with their floors",
    "beyond_controls": "Part 1's seven controls: could the leftover be an artefact?",
    "beyond_where": "Where does the map sit above or below what they predict?",
    "one_comparison": "What is a gene's rho with the map?",
    "one_comparison_detail": "One comparison step by step: the maps as measured, as "
    "ranks, for three genes",
    "spatial_null": "How large a rho do unrelated smooth maps give?",
    "spatial_null_detail": "The spatial null in detail: surrogates on a plane, two "
    "genes against their null",
    "top_genes": "Which genes follow the map, and what kind of maps are they?",
    "gene_ranking": "P9's 100 genes against the map, the Cacng8 - Gria1 gap in "
    "detail, the two genes with the tissue",
    "cacng8_gria1": "Does the map follow Cacng8 more closely than Gria1?",
    "between_within": "Do genes follow the map inside divisions, or only between them?",
    "between_within_detail": "Between and within divisions in detail: five genes "
    "division by division, the choice of null",
    "gene_kinds": "Do the genes that set surface receptor follow the map better than "
    "other genes?",
    "gene_sets": "The gene sets in detail: every set, the two contrasts named in "
    "advance, where each set comes from",
    "localisation": "The localisation test in detail: its null, the positive control, "
    "the matching, every test of the design",
    "leftover": "Does any gene follow what Gria1 and synapse density leave?",
    "leftover_genes": "Every gene and every gene set against the leftover",
    "ampa_family": "The AMPA receptor complex family against the leftover, member by "
    "member, and on the map",
    "autofluorescence": "Is the gene ranking the label's or the tissue's?",
    "autofluorescence_detail": "The autofluorescence map in detail: the genes' rho with "
    "each map, the counts past each null",
    "robustness": "Does the ranking depend on the choices made?",
    "robustness_detail": "The choices made, with every gene under four of them",
    "synaptome": "Does a measured synapse density change part 1?",
    "synaptome_detail": "The measured synapse density in detail: coverage, the two "
    "hemispheres, agreement with every map",
    "density_markers": "Which genes stand for synapse density, chosen without the map?",
    "green_channel": "Does the green channel report the tagged receptor, or the tissue?",
    "green_channel_detail": "The three channels in detail: one adult's raw planes, each "
    "channel against Gria1",
    "april_headline": "What is left of April's headline?",
    "april_headline_detail": "April's ten groups then and now, and against the null",
}


# the main figures of each part of the walk, in order
PARTS = (
    ("The question", ("overview",)),
    ("The inputs", ("structures", "genes")),
    (
        "Part 1: the map is not fully explained by Gria1 expression and synapse density",
        ("beyond", "beyond_where"),
    ),
    (
        "Part 2: what else it is: the genes that follow the map and its leftover",
        (
            "one_comparison",
            "spatial_null",
            "top_genes",
            "cacng8_gria1",
            "between_within",
            "gene_kinds",
            "leftover",
        ),
    ),
    (
        f"Controls ({FIGURE_NUMBERS['autofluorescence']} and "
        f"{FIGURE_NUMBERS['robustness']} of part 2, {FIGURE_NUMBERS['synaptome']} of "
        "part 1)",
        ("autofluorescence", "robustness", "synaptome"),
    ),
    ("The limit", ("green_channel",)),
    ("April's headline", ("april_headline",)),
)

# the detailed versions of each main figure, and what each adds
SUPPLEMENTS = {
    "genes": (
        (
            "genes_detail",
            "the section QC of P9's experiments, reliability against expression, and "
            "what was left out and repaired",
        ),
    ),
    "beyond": (
        (
            "beyond_budget",
            "the map against each term, every model's share, the four parts with "
            "their intervals and the weights, every draw of the calibration, and the "
            "main model under other folds and every check row, each with its own floor",
        ),
        (
            "beyond_controls",
            "the seven controls, one panel each: a gradient, structure size, single "
            "animals, naive against RWS, curvature, the whole gene table, the reading",
        ),
    ),
    "one_comparison": (
        (
            "one_comparison_detail",
            "the maps as measured and as ranks for nano, Cacng8, Gria1 and Aqp4, and "
            "the steps of one comparison",
        ),
    ),
    "spatial_null": (
        (
            "spatial_null_detail",
            "the nano map and three surrogates on a plane, and Cacng8's and Gria1's rho "
            "against their nulls",
        ),
    ),
    "top_genes": (
        (
            "gene_ranking",
            "P9's 100 genes one by one with their null bands and autofluorescence's "
            "rho, and the Cacng8 - Gria1 gap with the adults' interval and each pairing "
            "of Allen experiments",
        ),
    ),
    "between_within": (
        (
            "between_within_detail",
            "five genes division by division, the choice of null for the within rho, "
            "and the genes highest inside divisions",
        ),
    ),
    "gene_kinds": (
        (
            "gene_sets",
            "every set with its genes named, the two contrasts named in advance, and "
            "where each set comes from",
        ),
        (
            "localisation",
            "the label null, the positive control with both control pools, the "
            "matching on expression, and every test of the design",
        ),
    ),
    "leftover": (
        (
            "leftover_genes",
            "the genes closest to the leftover with their null bands, and every gene "
            "set against it",
        ),
        (
            "ampa_family",
            "each member of the family against the leftover, and the family on the "
            "map itself",
        ),
    ),
    "autofluorescence": (
        (
            "autofluorescence_detail",
            "the genes' rho with each map, how many pass each null at three "
            "thresholds, and the genes that pass",
        ),
    ),
    "robustness": (("robustness_detail", "every gene under four of the choices"),),
    "synaptome": (
        (
            "synaptome_detail",
            "the two hemispheres of the one mouse, and each density's agreement with "
            "the density term, Gria1 and the maps",
        ),
        (
            "density_markers",
            "the postsynaptic genes of the density term, chosen by their agreement with "
            "PSD95 puncta without the map, the AMPA-linked genes left out, and the "
            "choice on random halves of the structures",
        ),
    ),
    "green_channel": (
        (
            "green_channel_detail",
            "one adult's raw channels on a plane, and each channel against Gria1, with "
            "what is left of SEP once autofluorescence is out",
        ),
    ),
    "april_headline": (
        ("april_headline_detail", "April's ten violins beside today's, gene by gene"),
    ),
}

# the run script that draws each guided figure
DRAWN_BY = {
    "overview": "run_ish_overview.py",
    "structures": "run_structure_set.py",
    "genes": "run_ish_gene_table.py",
    "genes_detail": "run_ish_gene_table.py",
    "beyond": "run_beyond_figures.py",
    "beyond_budget": "run_beyond_figures.py",
    "beyond_controls": "run_beyond_figures.py",
    "beyond_where": "run_beyond_figures.py",
    "one_comparison": "run_ish_gene_ranking.py",
    "one_comparison_detail": "run_ish_gene_ranking.py",
    "spatial_null": "run_ish_gene_ranking.py",
    "spatial_null_detail": "run_ish_gene_ranking.py",
    "top_genes": "run_ish_top_genes.py",
    "gene_ranking": "run_ish_gene_ranking.py",
    "cacng8_gria1": "run_ish_top_genes.py",
    "between_within": "run_ish_divisions.py",
    "between_within_detail": "run_ish_divisions.py",
    "gene_kinds": "run_ish_gene_sets.py",
    "gene_sets": "run_ish_gene_sets.py",
    "localisation": "run_ish_gene_sets.py",
    "leftover": "run_ish_top_genes.py",
    "leftover_genes": "run_beyond_figures.py",
    "ampa_family": "run_ish_top_genes.py",
    "autofluorescence": "run_ish_gene_ranking.py",
    "autofluorescence_detail": "run_ish_gene_ranking.py",
    "robustness": "run_ish_robustness.py",
    "robustness_detail": "run_ish_robustness.py",
    "synaptome": "run_beyond_figures.py",
    "synaptome_detail": "run_synaptome.py",
    "density_markers": "run_density_markers.py",
    "green_channel": "run_sep_channel_check.py",
    "green_channel_detail": "run_sep_channel_check.py",
    "april_headline": "run_ish_overview.py",
    "april_headline_detail": "run_ish_overview.py",
}


def figure_file(key: str) -> str:
    """The file name of a guided figure: its number and its key, 07s_gene_ranking.png."""
    return f"{FIGURE_NUMBERS[key]}_{key}.png"


def figure_path(key: str) -> Path:
    """Where a guided figure is written: figures/ of the ISH line, its file name."""
    return FIGURES / figure_file(key)


def figure_ref(key: str) -> str:
    """How a figure names another: 'figure 06', 'figure 10s2'."""
    return f"figure {FIGURE_NUMBERS[key]}"
