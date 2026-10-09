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
                               documentation; figures 02, 02s
    16. run_ish_spatial_null   A7: surrogate maps and their checks
    17. run_ish_gene_ranking   analysis 1: each gene against the
                               map, the null, autofluorescence (A8),
                               the Cacng8 - Gria1 gap; figures 05,
                               06, 07s and 12, with their s
    18. run_ish_robustness     A3: the ranking under other choices;
                               figures 13, 13s
    19. run_ish_divisions      analysis 2 (A6): between or within
                               divisions; figures 09, 09s, gene sheets
    20. run_ish_gene_sets      analysis 3: kinds of genes; localisation
                               against matched controls; figures 10,
                               10s1, 10s2
    21. run_synaptome          the measured synapse density (network,
                               once): PSD95 puncta per structure,
                               its coverage of the fit; figure 14s1
    22. run_density_markers    the synapse-density genes, chosen by
                               PSD95 without the map (network,
                               once); figure 14s2
    23. run_beyond_density     analysis 4: what Gria1 and synapse
                               density leave; the leftover
    24. run_beyond_controls    seven attempts to break it
    25. run_beyond_calibration the same model on maps whose answer
                               is known; the check rows, each
                               with its own floor
    26. run_beyond_regression  the regression, per structure
    27. run_beyond_figures     figures 03, 03s1, 03s2, 04, 11s1, 14
    28. run_ish_top_genes      the genes that follow the map most,
                               characterised; Cacng8 and the AMPA
                               receptor complex family against the
                               leftover; figures 07, 08, 11, 11s2,
                               sheets
    29. run_sep_channel_check  analysis 5: what the green channel
                               reports; figures 15, 15s
    30. run_ish_overview       figures 00, 16 and 16s, the figure
                               index; the numbers for the text

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
    figures/02s_genes_detail.png    the same with the section QC and the repair

    python run_ish_gene_table.py [--offline]

--offline stops instead of asking mygene.info or the Gene Ontology for what the
caches lack.
"""

import argparse

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

from sepmap import config, plotting, structures
from sepmap.ish import gene_table, section_qc
from sepmap.ish import plotting as ish_plotting
from sepmap.ish.figure_index import figure_file, figure_path


def main(offline):
    """Print the settings, then measure, label and document every gene."""
    config.print_settings({"offline": offline})

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
    region = gene_table.region_table(experiments, flags)
    region.to_csv(gene_table.REGION_TABLE, index=False, float_format="%.6g")
    print(
        f"region means: {len(region):,} rows, {region['experiment_id'].nunique()} "
        f"experiments"
    )

    # experiments that measure too few structures, then the reliability and the
    # merged profile of each gene from the rest, at full precision
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
    numbers = gene_table.numbers_table(experiments, rel, members, labels, release)
    numbers.to_csv(gene_table.NUMBERS, index=False)
    print(f"gene table: {len(table)} experiments, {len(documentation)} genes documented")

    # figure 02 and its detailed version
    usable_genes = labels[labels["symbol"].isin(set(profiles["symbol"]))]
    sections = section_qc.load_section_qc()
    fig = ish_plotting.plot_genes(usable_genes, table, rel, save=figure_path("genes"))
    plt.close(fig)
    fig = ish_plotting.plot_genes_detail(
        usable_genes,
        table,
        sections,
        summary,
        rel,
        set(members["subunits"]),
        save=figure_path("genes_detail"),
    )
    plt.close(fig)
    drawn = [figure_file(k) for k in ("genes", "genes_detail")]
    print(f"figures: {', '.join(drawn)} in {structures.FIGURES}")


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
