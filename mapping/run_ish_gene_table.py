"""The gene table of the ISH line: region means after QC, profiles, labels, sets.

Python route, in run order (tools\\venv_atlas; run_closeup in tools\\venv_flat):
     1. run_per_mouse          per brain: tissue mask, backgrounds
     2. run_to_ccf             every brain on the adult CCF grid
     3. run_cohort             cohort mean, SD and n per voxel
     4. run_compare            young against adult: maps, table
        run_replot             redraws those maps (not in a full run)
     5. run_region_plot        region statistics, no warping
     6. run_region_groups      the same by system and by layer
     7. run_video              cohort videos
     8. run_video_compare      young beside adult, plane by plane
     9. run_closeup            close-ups and flatmaps (venv_flat)
    10. run_diagnostics        sheets that audit each step
    11. run_panel_build        390-gene ontology panel (network, cached)
    12. run_panel_fetch        its ISH grids (network, once)
    13. run_structure_set      A1: the declared structures, the
                               adult profiles; figure 01
    14. run_ish_section_qc     A2: the experiments of both panels
                               and the repair (network, once);
                               section QC; QC sheets
    15. run_ish_gene_table     A9: region means, the gene table,    <- this script
                               merged profiles, gene sets,
                               documentation; figure 02
    16. run_ish_spatial_null   A7: surrogate maps and their checks
    17. run_ish_gene_ranking   analysis 1: each gene against the
                               map, the null, autofluorescence (A8),
                               the Cacng8 - Gria1 gap; figures 05 to 07 and 12
    18. run_ish_robustness     A3: the ranking under other choices;
                               figure 13
    19. run_ish_divisions      analysis 2 (A6): between or within
                               divisions; figure 10, gene sheets
    20. run_ish_gene_sets      analysis 3: kinds of genes; localisation
                               against matched controls; figures 08, 09
    21. run_beyond_density     analysis 4: what receptor mRNA and
                               synaptic density leave; the leftover
    22. run_beyond_controls    seven attempts to break it
    23. run_beyond_calibration the same model on maps whose answer
                               is known
    24. run_beyond_regression  the regression, per structure
    25. run_beyond_figures     figures 03, 04 and 11
    26. run_sep_channel_check  analysis 5: what the green channel
                               reports; figure 14
    27. run_ish_overview       figures 00 and 15, the figure index;
                               the numbers for the text

Takes the experiments and the section flags of run_ish_section_qc, sets the flagged
sections missing, measures every experiment per structure, and gathers per gene its
reliability, its merged profile, its labels (P9's category, the ontology panel's
role, the name and GO annotation of mygene.info) and the gene sets of analysis 3
(the method is in sepmap/ish/gene_table.py and sepmap/ish/gene_sets.py). Writes, in
adult_v2/ish_analysis/ under the data root:

    tables/gene_region_table.csv    per experiment and structure: full and eroded
                                    mean energy, voxels with data, coverage
    tables/gene_table.csv           per experiment: panels, QC, use or the reason
                                    not, and its gene's labels and reliability
    tables/gene_profiles.csv        per gene and structure, the mean rank of its
                                    experiments, the profile every analysis reads
    tables/gene_documentation.csv   per gene, for the supplement (UTF-8 with BOM)
    tables/numbers_gene_table.csv   the numbers of this step, for the text
    cache/                          mygene records and go-basic.obo
    figures/02_genes.png            the genes and how good their maps are

    python run_ish_gene_table.py [--offline]

--offline stops instead of asking mygene.info or the Gene Ontology for what the
caches lack.
"""

import argparse

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from sepmap import config, plotting
from sepmap.ish import gene_table, regions, section_qc
from sepmap.ish import plotting as ish_plotting
from sepmap.volumes.per_mouse import structure_terms

OUT = config.DATA / "adult_v2" / "ish_analysis"


def numbers_table(experiments, rel, members, labels, release):
    """The numbers of this step that the text quotes, one row each."""
    used = experiments[~experiments["excluded"]]
    values = rel["reliability"].dropna()
    by_gene = rel.set_index("symbol")["reliability"]
    p9_genes = set(labels.loc[labels["p9_gene"], "symbol"])
    rows = [
        ("genes_listed", experiments["symbol"].nunique(), "genes of both panels"),
        ("experiments_listed", len(experiments), "experiments, repair included"),
        ("experiments_used", len(used), "experiments with a usable grid"),
        ("genes_with_profile", rel["symbol"].nunique(), "genes with a profile"),
        (
            "p9_genes_with_profile",
            len(p9_genes & set(rel["symbol"])),
            "of P9's 100 genes",
        ),
        ("genes_reliability", len(values), "genes measured more than once"),
        ("reliability_median", round(float(values.median()), 3), "median reliability"),
        ("reliability_q1", round(float(values.quantile(0.25)), 3), "its first quartile"),
        ("reliability_q3", round(float(values.quantile(0.75)), 3), "its third quartile"),
        ("reliability_below_0.3", int((values < 0.3).sum()), "genes below 0.3"),
        ("go_release", release, "go-basic.obo release of the ancestry"),
    ]
    for gene in ("Gria1", "Cacng8", "Dlg2"):
        value = by_gene.get(gene, np.nan)
        rows.append(
            (f"reliability_{gene}", round(float(value), 3), f"{gene}'s reliability")
        )
    for name, symbols in members.items():
        rows.append((f"set_{name}", len(symbols), f"genes in the set {name}"))
    for reason, n in (
        experiments.loc[experiments["excluded"], "exclude_reason"].value_counts().items()
    ):
        rows.append(("excluded", int(n), reason))
    return pd.DataFrame(rows, columns=["name", "value", "what"], dtype=object)


def main(offline):
    """Print the settings, then measure, label and document every gene."""
    config.print_settings({"offline": offline})
    tables = OUT / "tables"

    # the experiments, which can be used, and the sections to set missing
    experiments = gene_table.exclusion_reasons(gene_table.load_experiments())
    summary = section_qc.load_experiment_qc()
    flags = section_qc.flagged_sections(summary)
    n_used = int((~experiments["excluded"]).sum())
    print(
        f"experiments: {len(experiments)} listed, {n_used} usable, "
        f"{sum(len(v) for v in flags.values())} sections to set missing in "
        f"{len(flags)} of them"
    )

    # region means per experiment, flagged sections missing
    ann = regions.annotation_200()
    names, _, _ = structure_terms()
    eroded = regions.eroded_annotation(ann)
    region = pd.DataFrame(gene_table.region_rows(experiments, flags, ann, names, eroded))
    region.to_csv(gene_table.REGION_TABLE, index=False, float_format="%.6g")
    print(
        f"region means: {len(region):,} rows, {region['experiment_id'].nunique()} "
        f"experiments"
    )

    # experiments that measure too few structures, then the reliability and the
    # merged profile of each gene from the rest
    n_before = int(experiments["excluded"].sum())
    experiments = gene_table.coverage_exclusions(experiments, region)
    used = set(experiments.loc[~experiments["excluded"], "experiment_id"])
    n_few = int(experiments["excluded"].sum()) - n_before
    print(f"coverage: {n_few} experiments measure too few structures; {len(used)} used")
    per = gene_table.experiment_profiles(region[region["experiment_id"].isin(used)])
    rel = gene_table.gene_reliability(per)
    profiles = pd.DataFrame(gene_table.profile_rows(per))
    profiles.to_csv(gene_table.PROFILES, index=False, float_format="%.6g")
    n_rel = int(rel["reliability"].notna().sum())
    print(
        f"profiles: {profiles['symbol'].nunique()} genes; reliability for {n_rel}, "
        f"median {rel['reliability'].median():.3f}"
    )

    # labels: mygene records, GO with its ancestry, the gene sets
    genes = sorted(experiments["symbol"].unique())
    records = gene_table.mygene_records(genes, offline)
    parents, _, release = gene_table.load_obo(offline)
    labels, members = gene_table.gene_labels(experiments, records, parents)
    sizes = ", ".join(f"{k} {len(v)}" for k, v in members.items())
    n_records = sum(1 for r in records.values() if r)
    print(f"labels: {n_records} mygene records, GO {release}; sets: {sizes}")

    # the gene table, the documentation, the numbers
    table = gene_table.experiment_table(experiments, summary, region, labels, rel)
    table.to_csv(gene_table.GENE_TABLE, index=False)
    documentation = gene_table.documentation_table(labels, table, rel)
    documentation.to_csv(gene_table.DOCUMENTATION, index=False, encoding="utf-8-sig")
    numbers = numbers_table(experiments, rel, members, labels, release)
    numbers.to_csv(tables / "numbers_gene_table.csv", index=False)
    print(f"gene table: {len(table)} experiments, {len(documentation)} genes documented")

    # figure 02
    usable_genes = labels[labels["symbol"].isin(set(profiles["symbol"]))]
    sections = section_qc.load_section_qc()
    path = OUT / "figures" / ish_plotting.figure_file("genes")
    fig = ish_plotting.plot_genes(
        usable_genes,
        table,
        sections,
        summary,
        rel,
        set(members["subunits"]),
        save=path,
    )
    plt.close(fig)
    print(f"figure: {path}")


if __name__ == "__main__":
    # figures go to files, never to a window
    matplotlib.use("Agg")
    plotting.set_style()

    parser = argparse.ArgumentParser(description="the gene table of the ISH line (A9)")
    parser.add_argument(
        "--offline",
        action="store_true",
        help="never ask mygene.info or the Gene Ontology; stop instead",
    )
    args = parser.parse_args()
    main(args.offline)
