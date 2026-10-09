"""The check rows of part 1: the main model changed one thing at a time, with floors.

The main model of adult.beyond_density (Gria1 and synapse density, two straight
terms) gives the value quoted. Each check row of settings.toml [beyond.checks]
changes one thing of it and is reported beside it, never in its place, so a reader
sees how much each choice matters:

    curved            each term bent (x, x^2, x^3). A straight model leaves any
                      curvature of the map's relation to its terms in the leftover,
                      and the floor of a straight model holds none, so this row shows
                      how much of the main one is curvature. A check, not a bound:
                      terms bent further take a little more (control E)
    autofluorescence  + the autofluorescence of the same sections, the one predictor
                      measured in our own brains
    first_proposal    the density genes proposed before the rule (Dlg4, Homer1,
                      Camk2a; the rule excludes Dlg4 and Camk2a as AMPA-receptor-linked)
    marker_panel      the 11 synaptic marker genes of the marker panel ([beyond]
                      markers), presynaptic and postsynaptic
    psd_pc1           + the first component of the postsynaptic-density genes of the
                      ontology panel, in the density group: those measured on every
                      structure, and in the calibration those with two halves of
                      experiments, one list for both (the columns psd_genes)
    four_subunits     Gria1 to Gria4 in place of Gria1, each its own term
    psd95             the measured PSD95 punctum density as the density term, on the
                      structures it covers (adult.synaptome)
    psd95_main        the main model on those same structures, so the two density
                      measures meet on equal ground
    large             only the structures of at least beyond.large_voxels voxels of
                      the CCF annotation on the 200 um grid
    allen_grid        nano measured on the Allen 200 um grid with the structures
                      assigned as the ISH values are (adult.profiles), so that a
                      mismatch of grids cannot pass for a leftover

Each row runs on its own structures, the declared ones where its predictors are
measured (or the main model's, cut as the row says), with its ceiling recomputed
there; its shares are held out as the main model's. Each has its own floor: the
calibration of adult.beyond_calibration with the row's own model, on the row's
structures where both halves of the Allen experiments measure its genes, random
folds, and nano minus that floor with its interval from the same paired jackknife
over structures, and over spatial blocks. A floor misses the mismatch of what is
measured once, the same in both halves: the PSD95 density of the row psd95 (one
mouse), and any gene of a row with a single usable experiment (the column once names
them), so such a floor errs low by more than the main model's.

Writes, in adult_v2/ish_analysis/beyond/ under the data root:

    check_rows.csv              per check row: its structures, terms, ceiling, the
                                shares of abundance, density and the whole model,
                                the share left, the calibration's structures, nano
                                there, the floor, nano minus the floor with its 95%
                                intervals (over structures, over spatial blocks),
                                what is measured once, and the genes behind psd_pc1
    check_rows_calibration.csv  every row's calibration (as calibration.csv, random
                                folds, without the leftovers' replication)
    check_rows_jackknife.csv    every row's paired jackknife of nano and the floor
    check_rows_jackknife_blocks.csv  the same over spatial blocks

Run by run_beyond_calibration.py.
"""

import dataclasses

import numpy as np
import pandas as pd

from sepmap.adult import beyond_calibration, beyond_density, profiles, synaptome
from sepmap.config import SETTINGS

# the check rows by key, with their labels; the voxels of the row of large structures
BEYOND = SETTINGS["beyond"]
CHECKS = BEYOND["checks"]
ISH = SETTINGS["ish"]

CHECK_ROWS = beyond_density.OUT / "check_rows.csv"
CHECK_CALIBRATION = beyond_density.OUT / "check_rows_calibration.csv"
CHECK_JACKKNIFE = beyond_density.OUT / "check_rows_jackknife.csv"
CHECK_JACKKNIFE_BLOCKS = beyond_density.OUT / "check_rows_jackknife_blocks.csv"


# ===== The rows =====


def large_structures(structures: list[str]) -> list[str]:
    """Those of `structures` with beyond.large_voxels 200 um voxels in the CCF annotation.

    Counted on the annotation as ish.regions samples it onto the Allen grid, both
    hemispheres, as an ISH grid's structure is.
    """
    codes, label_names, _ = profiles.grid_codes()
    n = np.bincount(codes.ravel(), minlength=len(label_names))
    size = dict(zip(label_names, n))
    return [s for s in structures if size.get(s, 0) >= BEYOND["large_voxels"]]


def grid_inputs(inputs: beyond_density.Inputs) -> beyond_density.Inputs:
    """The inputs with nano as measured on the Allen 200 um grid (adult.profiles).

    On the inputs' structures where every adult has ish.min_voxels voxels of nano on
    the grid; the adults' zref there replaces nano, the rest is unchanged.
    """
    grid = profiles.load_grid_table()
    grid = grid[grid["n_voxels"] >= ISH["min_voxels"]]
    wide = grid.pivot(index="mouse", columns="structure", values="zref_nano")
    wide = wide.reindex(index=profiles.ADULTS)
    complete = set(wide.columns[wide.notna().all(axis=0)])
    kept = [s for s in inputs.structures if s in complete]
    sub = beyond_density.restrict(inputs, kept)
    return dataclasses.replace(sub, nano=wide[kept].to_numpy(float))


def row_inputs(
    terms: dict[str, tuple[str, ...]],
) -> beyond_density.Inputs:
    """The inputs on the declared structures where a model's genes are measured."""
    return beyond_density.load_inputs(beyond_density.model_genes(terms))


def check_models(
    inputs: beyond_density.Inputs,
) -> list[tuple[str, beyond_density.Inputs, dict[str, tuple[str, ...]], bool]]:
    """The check rows as (key, inputs, terms, bend), in the order of [beyond.checks].

    `inputs` are the main model's.
    """
    main = inputs.terms
    psd95 = (synaptome.MEASURE,)
    first = beyond_density.model_terms(("first_proposal",))
    markers = beyond_density.model_terms(("markers",))
    subunits = beyond_density.model_terms(abundance=beyond_density.SUBUNITS)
    measured = beyond_density.model_terms(psd95)
    models = {
        "curved": (inputs, main, True),
        "autofluorescence": (
            inputs,
            beyond_density.model_terms(autofluorescence=True),
            False,
        ),
        "first_proposal": (row_inputs(first), first, False),
        "marker_panel": (row_inputs(markers), markers, False),
        "psd_pc1": (inputs, beyond_density.model_terms(("density", "psd_pc1")), False),
        "four_subunits": (row_inputs(subunits), subunits, False),
        "psd95": (
            beyond_density.on_measured(row_inputs(measured), psd95),
            measured,
            False,
        ),
        "psd95_main": (beyond_density.on_measured(inputs, psd95), main, False),
        "large": (
            beyond_density.restrict(inputs, large_structures(inputs.structures)),
            main,
            False,
        ),
        "allen_grid": (grid_inputs(inputs), main, False),
    }
    return [(key, *models[key]) for key in CHECKS]


def terms_text(terms: dict[str, tuple[str, ...]], bend: bool) -> str:
    """A model's predictors in one line: 'Gria1; density, psd_pc1; autofluo; bent'."""
    groups = [", ".join(terms[g]) for g in beyond_density.ORDER if g in terms]
    return "; ".join(groups) + ("; bent" if bend else "")


def psd_gene_counts(
    inputs: beyond_density.Inputs,
    terms: dict[str, tuple[str, ...]],
    halves: dict,
    psd: list[str],
) -> tuple[float, float]:
    """The genes behind psd_pc1 in the row's fit and in its calibration; NaN without it.

    `psd` is the list the fit used (beyond_density.covariates_for).
    """
    if "psd_pc1" not in terms["density"]:
        return float("nan"), float("nan")
    cal = beyond_calibration.calibration_inputs(inputs, halves, terms)
    shared = beyond_calibration.shared_psd_genes(cal, halves)
    return float(len(psd)), float(len(shared))


def check_row(
    key: str,
    inputs: beyond_density.Inputs,
    model: tuple[dict[str, tuple[str, ...]], bool],
    halves: dict,
    split: list[str],
) -> tuple[dict, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """One check row: its shares, its floor and nano minus the floor with its intervals.

    `model` is the row's terms and whether they bend; `halves` and `split` are those
    of beyond_calibration.half_profiles. Returns the row, its calibration and its
    paired jackknives over structures and over spatial blocks, each with the row's
    key.
    """
    terms, bend = model
    splits = beyond_density.half_splits()
    _, explainable = beyond_density.ceiling(inputs.nano, splits)
    y = beyond_density.full_map(inputs.nano)
    cov, psd, _ = beyond_density.covariates_for(inputs)
    shares = beyond_density.model_shares(y, cov, terms, explainable, bend=bend)

    # the row's own floor, and nano minus it on the same subsamples, at random and
    # one spatial block at a time
    calibration = beyond_calibration.calibration_rows(
        inputs, halves, terms, bend, blocks=False, replication=False
    )
    jack = beyond_calibration.paired_jackknife(inputs, halves, terms, bend)
    blocks = beyond_calibration.paired_jackknife(inputs, halves, terms, bend, blocks=True)
    floor = beyond_calibration.left_summary(calibration)
    nano = beyond_calibration.nano_left(calibration)
    point = nano - floor["left_median"]
    lo, hi = beyond_calibration.jackknife_interval(jack, jack["nano_minus_floor"], point)
    block_lo, block_hi = beyond_calibration.jackknife_interval(
        blocks, blocks["nano_minus_floor"], point
    )
    n_psd, n_psd_cal = psd_gene_counts(inputs, terms, halves, psd)
    row = dict(
        key=key,
        check=CHECKS[key],
        n_structures=len(inputs.structures),
        terms=terms_text(terms, bend),
        ceiling=explainable,
        abundance=shares["abundance"],
        density_alone=shares["density"],
        model=shares["model"],
        left=shares["left"],
        cal_structures=int(calibration["n_structures"].iloc[0]),
        nano_cal_left=nano,
        floor=floor["left_median"],
        floor_lo=floor["left_lo"],
        floor_hi=floor["left_hi"],
        nano_minus_floor=point,
        nano_minus_floor_lo=lo,
        nano_minus_floor_hi=hi,
        nano_minus_floor_blocks_lo=block_lo,
        nano_minus_floor_blocks_hi=block_hi,
        once=" ".join(beyond_calibration.measured_once(split, terms)),
        psd_genes=n_psd,
        psd_genes_calibration=n_psd_cal,
    )
    for table in (calibration, jack, blocks):
        table.insert(0, "key", key)
    return row, calibration, jack, blocks


# ===== Reading back =====


def load_check_rows() -> pd.DataFrame:
    """check_rows.csv as run_beyond_calibration wrote it, empty text where blank."""
    if not CHECK_ROWS.exists():
        raise FileNotFoundError(
            f"{CHECK_ROWS} not found: run run_beyond_calibration.py first"
        )
    table = pd.read_csv(CHECK_ROWS, keep_default_na=False, na_values=[""])
    table["once"] = table["once"].fillna("")
    return table


def main() -> None:
    """Run every check row with its floor and write the four tables."""
    inputs = beyond_density.load_inputs()
    halves, split = beyond_calibration.half_profiles()
    rows, calibrations, jackknives, in_blocks = [], [], [], []
    print(
        f"\nthe check rows, beside the main model on {len(inputs.structures)} structures"
        " (95% over structures; over spatial blocks):"
    )
    for key, mine, terms, bend in check_models(inputs):
        row, calibration, jack, blocks = check_row(
            key, mine, (terms, bend), halves, split
        )
        rows.append(row)
        calibrations.append(calibration)
        jackknives.append(jack)
        in_blocks.append(blocks)
        once = f"; the same in both halves: {row['once']}" if row["once"] else ""
        print(
            f"  {key:17s} {row['n_structures']:4d} structures, left {row['left']:6.1%}; "
            f"on {row['cal_structures']} nano {row['nano_cal_left']:6.1%}, floor "
            f"{row['floor']:6.1%}, nano minus floor {row['nano_minus_floor']:+6.1%} "
            f"({row['nano_minus_floor_lo']:+.1%} to {row['nano_minus_floor_hi']:+.1%}; "
            f"{row['nano_minus_floor_blocks_lo']:+.1%} to "
            f"{row['nano_minus_floor_blocks_hi']:+.1%}){once}",
            flush=True,
        )
    table = pd.DataFrame(rows)
    table.to_csv(CHECK_ROWS, index=False)
    pd.concat(calibrations, ignore_index=True).to_csv(CHECK_CALIBRATION, index=False)
    pd.concat(jackknives, ignore_index=True).to_csv(CHECK_JACKKNIFE, index=False)
    pd.concat(in_blocks, ignore_index=True).to_csv(CHECK_JACKKNIFE_BLOCKS, index=False)
    print(
        f"  -> {CHECK_ROWS.name}, {CHECK_CALIBRATION.name}, {CHECK_JACKKNIFE.name}, "
        f"{CHECK_JACKKNIFE_BLOCKS.name}"
    )
