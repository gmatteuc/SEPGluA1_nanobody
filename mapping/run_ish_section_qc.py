"""Section quality of every ISH experiment of the gene table: flags and QC sheets.

Python route, in run order (tools\\venv_atlas; run_closeup and run_adult_layers in
tools\\venv_flat):
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
    14. run_ish_section_qc     A2: the experiments of both panels   <- this script
                               and the repair (network, once);
                               section QC; QC sheets
    15. run_ish_gene_table     A9: region means, the gene table,
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
                               its coverage of the fit; figure 03s4
    22. run_density_markers    the synapse-density genes, chosen by
                               PSD95 without the map (network,
                               once); figure 03s3
    23. run_beyond_density     analysis 4: what Gria1 and synapse
                               density leave; the leftover
    24. run_beyond_controls    seven attempts to break it
    25. run_beyond_calibration the same model on maps whose answer
                               is known; the check rows, each
                               with its own floor
    26. run_beyond_regression  the regression, per structure
    27. run_beyond_figures     figures 03, 03s1, 03s2, 04, 11s1
    28. run_ish_top_genes      the genes that follow the map most,
                               characterised; Cacng8 and the AMPA
                               receptor complex family against the
                               leftover; figures 07, 08, 11, 11s2,
                               sheets
    29. run_sep_channel_check  analysis 5: what the green channel
                               reports; figures 14, 14s
    30. run_ish_overview       figures 00, 15 and 15s, the figure
                               index; the numbers for the text
    31. run_adult_layers       the adult map by depth (venv_flat)

Lists every experiment of the two gene panels, adds the other Allen experiments of
the P9 genes whose own grid cannot be used (fetched once, cached), and judges every
section of every grid against its neighbours: a section five-fold dimmer than them
is flagged, to be set missing by the next step, unless the exceptions list
(mapping/ish_section_exceptions.csv) keeps it as a true absence. The methods are in
sepmap/ish/section_qc.py and sepmap/ish/gene_table.py. Writes, in
adult_v2/ish_analysis/ under the data root:

    tables/experiments.csv            every experiment, and which panel lists it
    tables/section_qc.csv             per experiment and section: voxels, median,
                                      local reference, status and reason
    tables/experiment_qc.csv          per experiment: sections judged, flagged,
                                      kept as absence
    tables/numbers_section_qc.csv     the numbers of this step, for the text
    figures/qc/00_flagged.png         every experiment with a flag, to review the
                                      exceptions from
    figures/qc/<gene>_<experiment>.png  one QC sheet per experiment (--sheets)

    python run_ish_section_qc.py [--sheets] [--redraw] [--offline]

--sheets draws the QC sheet of every experiment whose sheet is not there yet (about
15 minutes for all of them), --redraw draws them all again, --offline stops instead
of asking the Allen API for an experiment list the cache lacks.
"""

import argparse

import matplotlib
import matplotlib.pyplot as plt

from sepmap import config, plotting, structures
from sepmap.ish import gene_table, planes, regions, section_qc
from sepmap.ish.plotting import inputs as input_figures

# the QC sheets and the overview of the flags
QC_FIGURES = structures.FIGURES / "qc"


def draw_sheets(sections, summary, shape, redraw):
    """One QC sheet per usable experiment, those already drawn only with `redraw`.

    Returns how many were drawn.
    """
    template, labels = planes.sheet_atlas()
    usable = summary[summary["grid"] == "ok"]
    n_drawn = 0
    for i, r in enumerate(usable.to_dict("records"), 1):
        path = QC_FIGURES / f"{r['symbol']}_{r['experiment_id']}.png"
        if path.exists() and not redraw:
            continue
        table = sections[sections["experiment_id"] == r["experiment_id"]]
        grid = section_qc.read_grid(r["experiment_id"], shape)
        fig = input_figures.plot_qc_sheet(table, r, grid, template, labels, save=path)
        plt.close(fig)
        n_drawn += 1
        if i % 50 == 0:
            print(f"  {i}/{len(usable)} sheets", flush=True)
    return n_drawn


def main(sheets, redraw, offline):
    """Print the settings, list the experiments, judge their sections, draw."""
    config.print_settings({"sheets": sheets, "redraw": redraw, "offline": offline})
    structures.TABLES.mkdir(parents=True, exist_ok=True)

    # every experiment of both panels, and the repair
    experiments = gene_table.experiment_list(offline)
    experiments.to_csv(gene_table.EXPERIMENTS, index=False)
    print(
        f"experiments: {len(experiments)} of {experiments['symbol'].nunique()} genes "
        f"(P9 {int(experiments['p9_experiment'].sum())}, ontology "
        f"{int(experiments['ontology_experiment'].sum())}, repair "
        f"{int(experiments['repair_experiment'].sum())})"
    )

    # the reviewed exceptions, and the CCF brain on the 200 um grid
    exceptions = section_qc.load_exceptions()
    ann = regions.annotation_200()
    print(f"exceptions: {len(exceptions)} listed in {section_qc.EXCEPTIONS.name}")

    # every section of every grid
    sections, summary = section_qc.qc_tables(experiments, ann > 0, exceptions)
    sections.to_csv(section_qc.SECTION_QC, index=False)
    summary.to_csv(section_qc.EXPERIMENT_QC, index=False)
    ok = summary[summary["grid"] == "ok"]
    print(
        f"section QC: {len(ok)} grids judged ({len(summary) - len(ok)} not usable), "
        f"{int((ok['n_flagged'] > 0).sum())} with flags, {int(ok['n_flagged'].sum())} "
        f"sections flagged, {int(ok['n_absence_kept'].sum())} kept as absence"
    )
    numbers = section_qc.numbers_table(experiments, summary, sections)
    numbers.to_csv(section_qc.NUMBERS, index=False)

    # the overview of the flags, and with --sheets one sheet per experiment
    QC_FIGURES.mkdir(parents=True, exist_ok=True)
    p9_genes = set(experiments.loc[experiments["p9_experiment"], "symbol"])
    fig = input_figures.plot_flagged(
        sections, summary, p9_genes, save=QC_FIGURES / "00_flagged.png"
    )
    plt.close(fig)
    print(f"figure: {QC_FIGURES / '00_flagged.png'}")
    if sheets or redraw:
        n_drawn = draw_sheets(sections, summary, ann.shape, redraw)
        print(f"QC sheets: {n_drawn} drawn in {QC_FIGURES}")


if __name__ == "__main__":
    # figures go to files, never to a window
    matplotlib.use("Agg")
    plotting.set_style()

    parser = argparse.ArgumentParser(description="section QC of the ISH experiments (A2)")
    parser.add_argument(
        "--sheets", action="store_true", help="draw the missing QC sheets"
    )
    parser.add_argument("--redraw", action="store_true", help="draw every QC sheet again")
    parser.add_argument(
        "--offline", action="store_true", help="never ask the Allen API; stop instead"
    )
    args = parser.parse_args()
    main(args.sheets, args.redraw, args.offline)
