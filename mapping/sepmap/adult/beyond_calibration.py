"""The model of analysis 4 on maps whose answer is known: the calibration floor.

A leftover is what the model does not predict. Part of it can be real, and part is
the predictors measuring something slightly different from what the map is made of:
each Allen experiment is one P56 mouse on a 200 um grid, and two experiments of one
gene agree at a median rho of about 0.7. So even a map made of nothing but receptor
mRNA and synaptic density would leave something, because the Allen maps that
predict it are not the maps it is made of. The calibration measures how much, with
the production model of adult.beyond_density:

    two halves    each gene's Allen experiments (those the gene table uses) are
                  split in two, alternately by id (the 1st, 3rd, ... against the
                  2nd, 4th, ...), and each half is merged as the gene table merges
                  (ish.reliability.merge); a gene measured once is in both halves.
                  The predictors are built from each half (build_covariates);
                  autofluorescence, measured in our own brains, is the same in both
    known maps    abundance and density   the production model fitted to the nano
                                          map with one half's predictors: exactly
                                          what receptor mRNA and synaptic density
                                          predict, with nano's own weights and
                                          curvature
                  Gria1 mRNA              Gria1's profile from one half: a map of
                                          one gene's mRNA, and nothing else
    animals       ten made-up adults per map: the known map (ranked, scaled to unit
                  variance) plus noise independent per structure and adult, of
                  variance 5 (1 / h - 1), so that two halves of five agree as nano's
                  halves do (h, their mean agreement)
    analysis      the production model with the other half's predictors: the
                  ceiling, the CV R2, the share left, and how well the two
                  half-cohort leftovers agree; both ways round, with
                  beyond_calibration.n_noise draws of the animals each way
    nano          the real map with each half's predictors, read on the same footing
                  as the known maps, and with the merged predictors of the
                  production run (predictors_from "both") on the same structures
    structures    those of beyond_density where both halves have every subunit and
                  marker gene: a structure that only one half measures drops out

A known map that is exactly abundance and density gives the floor: the share any
map would leave from one Allen map disagreeing with another. The Gria1 map shows how
much a map of one gene's mRNA leaves when the model knows that gene only from its
other experiments; its leftover replicates across the made-up animals too, which is
why replication alone cannot tell biology from ISH mismatch.

The noise is independent between structures, where real animals deviate smoothly;
the floor calibrates a size, it does not model the cohort.

Writes, in adult_v2/ish_analysis/beyond/ under the data root:

    calibration.csv   one row per map, direction and draw (draw -1 for the real map):
                      the structures, the half agreement, the ceiling, the CV R2, the
                      share left and the leftover's replication

Run by run_beyond_calibration.py.
"""

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from sepmap.adult import beyond_density as bd
from sepmap.config import SETTINGS
from sepmap.ish import gene_table, reliability

# the draws of animal noise per known map and direction
BEYOND_CALIBRATION = SETTINGS["beyond_calibration"]

CALIBRATION = bd.OUT / "calibration.csv"

# the known maps, as calibration.csv names them
NANO = "nano"
ABUNDANCE_DENSITY = "abundance and density"
GRIA1 = "Gria1 mRNA"

# the halves, as calibration.csv names them, and the predictors of the production
# run, both halves merged
HALVES = ("A", "B")
MERGED = "both"


def experiment_halves(
    per_experiment: dict[str, dict[str, dict[str, float]]],
) -> tuple[dict, dict, list[str]]:
    """Each gene's merged profile from half A and from half B of its experiments.

    `per_experiment` is {gene: {experiment: {structure: energy}}}. Experiments are
    taken in order of id, alternately into A and B; a gene measured once has the
    same profile in both. Returns the two {gene: {structure: mean rank}} and the
    genes whose halves differ.
    """
    a, b, split = {}, {}, []
    for gene, experiments in sorted(per_experiment.items()):
        ids = sorted(experiments, key=int)
        first = {e: experiments[e] for e in ids[0::2]}
        second = {e: experiments[e] for e in ids[1::2]} or first
        a[gene] = reliability.merge(first)[0]
        b[gene] = reliability.merge(second)[0]
        if len(ids) > 1:
            split.append(gene)
    return a, b, split


def standardised(values: np.ndarray) -> np.ndarray:
    """Ranks scaled to mean 0 and variance 1."""
    r = rankdata(values)
    return (r - r.mean()) / r.std()


def fake_cohort(
    truth: np.ndarray, half_agreement: float, rng: np.random.Generator
) -> np.ndarray:
    """Ten made-up adults, adults x structures: the known map plus animal noise.

    The noise has variance 5 (1 / h - 1) per adult and structure, so that the mean
    of five adults agrees with the mean of five others at about h; `truth` is
    standardised first.
    """
    sd = np.sqrt(bd.BEYOND["half"] * (1 / half_agreement - 1))
    z = standardised(truth)
    return z + rng.normal(0, sd, size=(len(bd.ADULTS), len(z)))


def analyse(
    cohort: np.ndarray,
    xs: list[np.ndarray],
    splits: list[tuple[list[int], list[int]]],
) -> dict:
    """The analysis of beyond_density on one cohort (adults x structures).

    Returns the half agreement, the ceiling, the CV R2, the share left and the
    leftover's replication, the mean over every split.
    """
    agreement, explainable = bd.ceiling(cohort, splits)
    cv = bd.cv_r2(bd.full_map(cohort), xs)
    return dict(
        half_agreement=float(np.mean(agreement)),
        ceiling=explainable,
        cv_r2=cv,
        left=1 - cv / explainable,
        replication=float(np.mean(bd.leftover_agreement(cohort, xs, splits))),
    )


def calibration_structures(
    structures: list[str], halves: dict[str, dict[str, dict[str, float]]]
) -> list[str]:
    """The structures of the fit where both halves have every subunit and marker."""
    needed = bd.SUBUNITS + bd.MARKERS
    return [
        s for s in structures if all(s in halves[h][g] for h in HALVES for g in needed)
    ]


def half_covariates(
    inputs: bd.Inputs, halves: dict, structures: list[str]
) -> dict[str, dict[str, np.ndarray]]:
    """The predictors built from each half of the experiments, on `structures`."""
    columns = [inputs.structures.index(s) for s in structures]
    auto = inputs.auto[:, columns]
    return {
        h: bd.build_covariates(halves[h], inputs.role, auto, structures)[0]
        for h in HALVES
    }


def calibration_rows(
    inputs: bd.Inputs,
    halves: dict[str, dict[str, dict[str, float]]],
    n_noise: int | None = None,
    seed: int = 0,
) -> pd.DataFrame:
    """calibration.csv: the real map and the known maps, analysed with each half.

    Each known map is built from one half and analysed with the other half's
    predictors, both ways round, n_noise draws each; every (map, direction) has a
    generator of its own, so a draw does not depend on which others ran.
    """
    if n_noise is None:
        n_noise = BEYOND_CALIBRATION["n_noise"]
    structures = calibration_structures(inputs.structures, halves)
    columns = [inputs.structures.index(s) for s in structures]
    nano = inputs.nano[:, columns]
    cov = half_covariates(inputs, halves, structures)
    splits = bd.half_splits()
    y = bd.full_map(nano)
    h_nano = float(np.mean(bd.ceiling(nano, splits)[0]))

    # the real map with each half's predictors, and with the merged ones of the
    # production run on the same structures
    merged, _, _ = bd.build_covariates(
        inputs.expr, inputs.role, inputs.auto[:, columns], structures
    )
    rows = []
    for h, covariates in [(h, cov[h]) for h in HALVES] + [(MERGED, merged)]:
        row = analyse(nano, bd.model(covariates), splits)
        rows.append(dict(map=NANO, truth_from="", predictors_from=h, draw=-1, **row))

    # the known maps: each from one half, analysed with the other's predictors
    known = {}
    for h in HALVES:
        xs = bd.model(cov[h])
        known[(ABUNDANCE_DENSITY, h)] = y - bd.residual(y, xs)
        known[(GRIA1, h)] = cov[h]["Gria1"]
    children = np.random.SeedSequence(seed).spawn(len(known))
    for ((kind, truth_from), truth), child in zip(known.items(), children):
        other = HALVES[1 - HALVES.index(truth_from)]
        xs = bd.model(cov[other])
        rng = np.random.default_rng(child)
        for draw in range(n_noise):
            cohort = fake_cohort(truth, h_nano, rng)
            row = analyse(cohort, xs, splits)
            rows.append(
                dict(
                    map=kind,
                    truth_from=truth_from,
                    predictors_from=other,
                    draw=draw,
                    **row,
                )
            )
    out = pd.DataFrame(rows)
    out.insert(1, "n_structures", len(structures))
    return out


def floor(table: pd.DataFrame, kind: str = ABUNDANCE_DENSITY) -> dict:
    """The median and range of the share left, and of the replication, for one map."""
    mine = table[table["map"] == kind]
    return dict(
        left_median=float(mine["left"].median()),
        left_lo=float(mine["left"].min()),
        left_hi=float(mine["left"].max()),
        replication_median=float(mine["replication"].median()),
        replication_lo=float(mine["replication"].min()),
        replication_hi=float(mine["replication"].max()),
    )


def load_calibration() -> pd.DataFrame:
    """The table run_beyond_calibration wrote."""
    if not CALIBRATION.exists():
        raise FileNotFoundError(
            f"{CALIBRATION} not found: run run_beyond_calibration.py first"
        )
    return pd.read_csv(CALIBRATION, keep_default_na=False, na_values=[""])


def per_experiment_profiles() -> dict[str, dict[str, dict[str, float]]]:
    """{gene: {experiment: {structure: energy}}} of the experiments the table uses."""
    table = gene_table.load_gene_table()
    region = gene_table.load_region_table()
    used = set(table.loc[~table["excluded"], "experiment_id"])
    return gene_table.experiment_profiles(region[region["experiment_id"].isin(used)])


def main() -> None:
    """Split the experiments, run the known maps and the real one, write the table."""
    bd.OUT.mkdir(parents=True, exist_ok=True)
    inputs = bd.load_inputs()
    a, b, split = experiment_halves(per_experiment_profiles())
    halves = {"A": a, "B": b}
    needed = bd.SUBUNITS + bd.MARKERS
    print(
        f"{len(split)} genes have two halves of experiments, "
        f"{sum(g in split for g in needed)} of the {len(needed)} subunit and marker "
        "genes; the others are the same in both"
    )
    table = calibration_rows(inputs, halves)
    table.to_csv(CALIBRATION, index=False)
    print(f"{int(table['n_structures'].iloc[0])} structures -> {CALIBRATION}")
    for r in table[table["map"] == NANO].itertuples():
        print(
            f"  nano, predictors from {r.predictors_from:4s}  left {r.left:6.1%}, "
            f"leftover replicates {r.replication:.3f}"
        )
    for kind in (ABUNDANCE_DENSITY, GRIA1):
        mine = table[table["map"] == kind]
        print(
            f"  {kind:22s} left {mine['left'].median():6.1%} "
            f"({mine['left'].min():.1%} to {mine['left'].max():.1%}), "
            f"leftover replicates {mine['replication'].median():.3f} "
            f"({mine['replication'].min():.3f} to {mine['replication'].max():.3f}), "
            f"{len(mine)} rows"
        )
