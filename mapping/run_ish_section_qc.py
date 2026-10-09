"""Section quality of every ISH experiment of the gene table: flags and QC sheets.

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
    14. run_ish_section_qc     A2: the experiments of both panels   <- this script
                               and the repair (network, once);
                               section QC; QC sheets
    15. run_ish_gene_table     A9: region means, the gene table,
                               merged profiles, gene sets,
                               documentation; figure 02
    16. run_ish_spatial_null   A7: surrogate maps and their checks
    17. run_ish_gene_ranking   analysis 1: each gene against the
                               map, the null, autofluorescence (A8),
                               the Cacng8 - Gria1 gap; figures 05 to 07 and 15
    18. run_ish_robustness     A3: the ranking under other choices;
                               figure 16
    19. run_ish_divisions      analysis 2 (A6): between or within
                               divisions; figure 10, gene sheets
    20. run_ish_gene_sets      analysis 3: kinds of genes; localisation
                               against matched controls; figures 08, 09
    21. run_synaptome          the measured synapse density (network,
                               once): PSD95 puncta per structure,
                               its coverage of the fit
    22. run_beyond_density     analysis 4: what Gria1 and synapse
                               density leave; the leftover
    23. run_beyond_controls    seven attempts to break it
    24. run_beyond_calibration the same model on maps whose answer
                               is known
    25. run_beyond_regression  the regression, per structure
    26. run_beyond_figures     figures 03, 04 and 11
    27. run_ish_top_genes      the genes that follow the map most,
                               characterised; Cacng8 and the AMPA
                               receptor complex family against the
                               leftover; figures 12 to 14, sheets
    28. run_sep_channel_check  analysis 5: what the green channel
                               reports; figure 17
    29. run_ish_overview       figures 00 and 18, the figure index;
                               the numbers for the text

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
import nibabel as nib
import numpy as np
import pandas as pd

from sepmap import config, plotting, structures
from sepmap.ish import gene_table, regions, section_qc
from sepmap.ish import plotting as ish_plotting
from sepmap.volumes.per_mouse import structure_terms

OUT = config.DATA / "adult_v2" / "ish_analysis"


def numbers_table(experiments, summary, sections):
    """The numbers of this step that the text quotes, one row each."""
    ok = summary[summary["grid"] == "ok"]
    reasons = sections["reason"].fillna("").astype(str)
    steps = sections[
        (sections["status"] == section_qc.OK)
        & reasons.str.startswith("dim against one side only")
    ]
    rows = [
        ("experiments_listed", len(experiments), "experiments of both panels and repair"),
        ("genes_listed", experiments["symbol"].nunique(), "genes of both panels"),
        (
            "experiments_repair",
            int(experiments["repair_experiment"].sum()),
            "experiments added by the repair",
        ),
        ("experiments_with_grid", len(ok), "experiments with a usable grid"),
        (
            "experiments_flagged",
            int((ok["n_flagged"] > 0).sum()),
            "experiments with at least one section set missing",
        ),
        ("sections_flagged", int(ok["n_flagged"].sum()), "sections set missing"),
        (
            "sections_absence_kept",
            int(ok["n_absence_kept"].sum()),
            "dim sections kept as true absence (exceptions list)",
        ),
        ("sections_judged", int(ok["n_judged"].sum()), "sections judged"),
        (
            "sections_faint",
            int((sections["status"] == section_qc.FAINT).sum()),
            "sections not judged: neighbours at the noise level",
        ),
        (
            "sections_step",
            len(steps),
            "dim against one side only: kept, a step in expression",
        ),
        (
            "experiments_step",
            steps["experiment_id"].nunique(),
            "experiments with a section kept at a step",
        ),
        (
            "experiments_near_zero",
            int(ok["near_zero"].astype(str).eq("True").sum()),
            "experiments whose median section is at the noise level",
        ),
    ]
    p9 = experiments.loc[experiments["p9_experiment"], "experiment_id"]
    p9_ok = ok[ok["experiment_id"].isin(set(p9))]
    rows.append(
        (
            "p9_experiments_flagged",
            int((p9_ok["n_flagged"] > 0).sum()),
            "P9's own experiments with a section set missing",
        )
    )
    return pd.DataFrame(rows, columns=["name", "value", "what"], dtype=object)


def draw_sheets(sections, summary, shape, redraw):
    """One QC sheet per experiment with a usable grid; those drawn only with `redraw`.

    The template and the structures are taken on the 20 um grid, ten times finer
    than the ISH grid, so the sheets draw the borders of structures.
    """
    folder = OUT / "figures" / "qc"
    template = nib.load(config.DATA / "atlas" / "average_template_10.nii.gz")
    template = np.asarray(template.dataobj)[::2, ::2, ::2]
    ann = np.asarray(nib.load(config.DATA / "atlas" / "annotation_10.nii.gz").dataobj)
    names, _, _ = structure_terms()
    labels, _ = structures.name_volume(ann[::2, ::2, ::2], names)
    usable = summary[summary["grid"] == "ok"]
    n_drawn = 0
    for i, r in enumerate(usable.to_dict("records"), 1):
        path = folder / f"{r['symbol']}_{r['experiment_id']}.png"
        if path.exists() and not redraw:
            continue
        table = sections[sections["experiment_id"] == r["experiment_id"]]
        grid = section_qc.read_grid(r["experiment_id"], shape)
        fig = ish_plotting.plot_qc_sheet(table, r, grid, template, labels, save=path)
        plt.close(fig)
        n_drawn += 1
        if i % 50 == 0:
            print(f"  {i}/{len(usable)} sheets", flush=True)
    return n_drawn


def main(sheets, redraw, offline):
    """Print the settings, list the experiments, judge their sections, draw."""
    config.print_settings({"sheets": sheets, "redraw": redraw, "offline": offline})
    tables = OUT / "tables"
    tables.mkdir(parents=True, exist_ok=True)

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
    numbers = numbers_table(experiments, summary, sections)
    numbers.to_csv(tables / "numbers_section_qc.csv", index=False)

    # the overview of the flags, and with --sheets one sheet per experiment
    folder = OUT / "figures" / "qc"
    folder.mkdir(parents=True, exist_ok=True)
    p9_genes = set(experiments.loc[experiments["p9_experiment"], "symbol"])
    fig = ish_plotting.plot_flagged(
        sections, summary, p9_genes, save=folder / "00_flagged.png"
    )
    plt.close(fig)
    print(f"figure: {folder / '00_flagged.png'}")
    if sheets or redraw:
        n_drawn = draw_sheets(sections, summary, ann.shape, redraw)
        print(f"QC sheets: {n_drawn} drawn in {folder}")


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
