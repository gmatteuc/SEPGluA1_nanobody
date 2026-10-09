"""The model of analysis 4 on maps whose answer is known: the calibration floor.

A leftover is what the model does not predict. Part of it can be real, and part is
the predictors measuring something slightly different from what the map is made of:
each Allen experiment is one P56 mouse on a 200 um grid, and two experiments of one
gene agree only in part (their reliability, ish.gene_table). So even a map made of
nothing but receptor mRNA and synapse density would leave something, because the
Allen maps that predict it are not the maps it is made of. The calibration measures
how much, with the main model of adult.beyond_density (Gria1, synapse density,
autofluorescence):

    two halves    each gene's Allen experiments (those the gene table uses) are
                  split in two, alternately by id (the 1st, 3rd, ... against the
                  2nd, 4th, ...), and each half is merged as the gene table merges
                  (ish.reliability.merge); a gene measured once is in both halves,
                  so its mismatch is not in the floor, which therefore errs low
                  (measured_once; numbers_beyond.csv counts them). The predictors
                  are built from each half (build_covariates); autofluorescence,
                  measured in our own brains, is the same in both, and so is the
                  measured PSD95 density when the main model uses it: one mouse,
                  measured once, like a marker with one experiment
    known maps    abundance and density   the main model fitted to the nano map
                                          with one half's predictors: exactly what
                                          Gria1 and synapse density predict, with
                                          nano's own weights and curvature
                  Gria1 mRNA              Gria1's profile from one half: a map of
                                          one gene's mRNA, and nothing else
    animals       ten made-up adults per map: the known map (ranked, scaled to unit
                  variance) plus noise independent per structure and adult, of
                  variance 5 (1 / h - 1), so that two halves of five agree as nano's
                  halves do (h, their mean agreement)
    analysis      the main model with the other half's predictors: the
                  ceiling, the CV R2, the share left, and how well the two
                  half-cohort leftovers agree; both ways round, with
                  beyond_calibration.n_noise draws of the animals each way
    nano          the real map with each half's predictors, read on the same footing
                  as the known maps, and with the merged predictors of the
                  production run (predictors_from "both") on the same structures
    structures    those of the main model where both halves have every subunit
                  and marker gene, as the fit's own rule asks of the merged
                  profiles: a structure that only one half measures drops out
    folds         every analysis twice: with the production folds (random, averaged
                  over beyond.cv_repeats shufflings) and with folds of spatial
                  blocks (beyond_density.block_labels), so the variant of blocked
                  folds has its own floor

A known map that is exactly Gria1 and synapse density gives the floor: the share
any map would leave from one Allen map disagreeing with another. The Gria1 map, the
benchmark, shows how much a map of one gene's mRNA leaves when the model knows that
gene only from its other experiments; its leftover replicates across the made-up
animals too, which is why replication alone cannot tell biology from ISH mismatch.

The noise is independent between structures, where real animals deviate smoothly;
the floor calibrates a size, it does not model the cohort. And it holds only the
mismatch of one Allen map with another, not that of Allen's P56 mice with these
brains (strain, age, the 200 um grid against 20 um masks, registration), so it
errs low there too.

The intervals of part 1 come from resampling structures, not animals: the claim is
about where in the brain the map departs from prediction, so what would differ in a
repeat is which structures were measurable. A delete-d jackknife: each of
beyond_calibration.n_jackknife subsamples leaves out a share
beyond_calibration.jackknife_share of the structures (jackknife_subsamples),
rebuilds every predictor on the structures kept as the production run builds them,
takes the ceiling from beyond_calibration.jackknife_splits splits of the adults, and
scores the model by the same cross-validation; the interval is the value plus or
minus 1.96 of its jackknife SD (jackknife_sd). A subsample draws no structure twice,
so no structure sits in a training and a test fold at once, as a bootstrap would let
it. Values compared with each other meet the same subsamples, so their difference
has an interval of its own: the nano map against the floor and the Gria1 map on the
calibration structures (paired_jackknife), and the check rows that change only the
density measure on the structures PSD95 covers (jackknife_density_rows), since the
interval of each alone is far wider than the differences between them.

Writes, in adult_v2/ish_analysis/beyond/ under the data root:

    calibration.csv   one row per folds, map, direction and draw (draw -1 for the
                      real map): the structures, the half agreement, the ceiling,
                      the CV R2, the share left and the leftover's replication
    calibration_jackknife.csv  per subsample of the jackknife: the nano map's share
                      left, the floor's, the Gria1 map's, and the differences

Run by run_beyond_calibration.py; adult.beyond_figures runs the jackknives of the
budget and of the density rows.
"""

from collections.abc import Iterator

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from sepmap.adult import beyond_density, synaptome
from sepmap.config import SETTINGS
from sepmap.ish import gene_table, reliability
from sepmap.structures import centroid_xyz

# the draws of animal noise per known map and direction; the jackknife's size; the
# gene of the stained protein, whose map is the benchmark
BEYOND_CALIBRATION = SETTINGS["beyond_calibration"]
ISH = SETTINGS["ish"]

CALIBRATION = beyond_density.OUT / "calibration.csv"
JACKKNIFE = beyond_density.OUT / "calibration_jackknife.csv"

# the maps, as calibration.csv names them: the real one, the floor (a map made of
# Gria1 and synapse density) and the benchmark (a map of one Gria1 experiment)
NANO = "nano"
FLOOR_MAP = "abundance and density"
GRIA1_MAP = "Gria1 mRNA"

# the halves, as calibration.csv names them, and the predictors of the production
# run, both halves merged
HALVES = ("A", "B")
MERGED = "both"

# the two kinds of folds, as calibration.csv names them
RANDOM = "random"
BLOCKS = "spatial blocks"

# the check rows whose density measures meet on the structures PSD95 covers, by their
# key in variants.csv, with their density terms
DENSITY_ROWS = {
    "psd95": beyond_density.MEASURED,
    "panel": beyond_density.PANEL,
    "psd95_and_panel": beyond_density.MEASURED + beyond_density.PANEL,
}

# the differences between those rows that the text quotes: (name, first, second,
# column), column "left" or "alone" (density alone)
DENSITY_CONTRASTS = (
    ("psd95_minus_panel_left", "psd95", "panel", "left"),
    ("both_minus_panel_left", "psd95_and_panel", "panel", "left"),
    ("psd95_minus_panel_alone", "psd95", "panel", "alone"),
)


# ===== The halves and the known maps =====


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


def half_profiles() -> tuple[dict[str, dict], list[str]]:
    """experiment_halves of the experiments the gene table uses, by half name.

    Returns {half: {gene: {structure: mean rank}}} and the genes whose halves
    differ.
    """
    table = gene_table.load_gene_table()
    a, b, split = experiment_halves(gene_table.used_experiment_profiles(table))
    return dict(zip(HALVES, (a, b))), split


def other_half(half: str) -> str:
    """The half of the experiments that is not `half`."""
    return HALVES[1 - HALVES.index(half)]


def standardised(values: np.ndarray) -> np.ndarray:
    """Ranks scaled to mean 0 and variance 1."""
    r = rankdata(values)
    return (r - r.mean()) / r.std()


def fake_cohort(
    truth: np.ndarray,
    half_agreement: float,
    rng: np.random.Generator | None = None,
    noise: np.ndarray | None = None,
) -> np.ndarray:
    """Ten made-up adults, adults x structures: the known map plus animal noise.

    The noise has variance 5 (1 / h - 1) per adult and structure, so that the mean
    of five adults agrees with the mean of five others at about h; `truth` is
    standardised first. The noise is drawn with `rng`, or scaled from `noise`
    (standard normal, adults x structures) when given.
    """
    sd = np.sqrt(beyond_density.BEYOND["half"] * (1 / half_agreement - 1))
    z = standardised(truth)
    if noise is None:
        noise = rng.standard_normal((len(beyond_density.ADULTS), len(z)))
    return z + sd * noise


def analyse(
    cohort: np.ndarray,
    columns: list[np.ndarray],
    splits: list[tuple[list[int], list[int]]],
    labels: list[np.ndarray] | None = None,
    replication: bool = True,
) -> dict:
    """The analysis of beyond_density on one cohort (adults x structures).

    Returns the half agreement, the ceiling, the CV R2 (folds `labels`, the
    production folds by default), the share left and, unless `replication` is
    False, the leftover's replication, the mean over every split.
    """
    agreement, explainable = beyond_density.ceiling(cohort, splits)
    cv = beyond_density.cv_r2(beyond_density.full_map(cohort), columns, labels)
    out = dict(
        half_agreement=float(np.mean(agreement)),
        ceiling=explainable,
        cv_r2=cv,
        left=1 - cv / explainable,
        replication=np.nan,
    )
    if replication:
        agree = beyond_density.leftover_agreement(cohort, columns, splits)
        out["replication"] = float(np.mean(agree))
    return out


def calibration_structures(
    structures: list[str], halves: dict[str, dict[str, dict[str, float]]]
) -> list[str]:
    """The structures of the fit where both halves have every subunit and marker."""
    needed = beyond_density.SUBUNITS + beyond_density.MARKERS
    return [
        s for s in structures if all(s in halves[h][g] for h in HALVES for g in needed)
    ]


def half_covariates(
    inputs: beyond_density.Inputs, halves: dict
) -> dict[str, dict[str, np.ndarray]]:
    """The predictors built from each half of the experiments, on the inputs' structures.

    Autofluorescence and the measured densities are the same in both.
    """
    return {
        h: beyond_density.build_covariates(
            halves[h],
            inputs.role,
            inputs.auto,
            inputs.structures,
            synapses=inputs.synapses,
        )[0]
        for h in HALVES
    }


def known_maps(
    y: np.ndarray,
    cov: dict[str, dict[str, np.ndarray]],
    terms: dict[str, tuple[str, ...]],
) -> dict:
    """The known maps, {(kind, half the truth comes from): ranks}.

    Abundance and density: the model `terms` fitted to the nano map `y` with one
    half's predictors, its fitted values; Gria1 mRNA: Gria1's ranks from one half.
    """
    out = {}
    for h in HALVES:
        columns = beyond_density.model_columns(cov[h], terms)
        out[(FLOOR_MAP, h)] = y - beyond_density.residual(y, columns)
        out[(GRIA1_MAP, h)] = cov[h][ISH["control_gene"]]
    return out


def calibration_rows(
    inputs: beyond_density.Inputs,
    halves: dict[str, dict[str, dict[str, float]]],
    n_noise: int | None = None,
    seed: int = 0,
) -> pd.DataFrame:
    """calibration.csv: the real map and the known maps, analysed with each half.

    The main model throughout (inputs.terms). Each known map is built from one
    half and analysed with the other half's predictors, both ways round, n_noise
    draws each, with random folds and with folds of spatial blocks; every (folds,
    map, direction) has a generator of its own, so a draw does not depend on which
    others ran.
    """
    if n_noise is None:
        n_noise = BEYOND_CALIBRATION["n_noise"]
    terms = inputs.terms
    cal = beyond_density.restrict(
        inputs, calibration_structures(inputs.structures, halves)
    )
    nano = cal.nano
    cov = half_covariates(cal, halves)
    splits = beyond_density.half_splits()
    y = beyond_density.full_map(nano)
    h_nano = float(np.mean(beyond_density.ceiling(nano, splits)[0]))
    folds = {
        RANDOM: None,
        BLOCKS: beyond_density.block_labels(centroid_xyz(cal.structures)),
    }

    # the real map with each half's predictors, and with the merged ones of the
    # production run on the same structures
    merged, _, _ = beyond_density.covariates_for(cal)
    rows = []
    known = known_maps(y, cov, terms)
    children = np.random.SeedSequence(seed).spawn(len(folds) * len(known))
    for f, (kind_of_folds, labels) in enumerate(folds.items()):
        for h, covariates in [(h, cov[h]) for h in HALVES] + [(MERGED, merged)]:
            columns = beyond_density.model_columns(covariates, terms)
            row = analyse(nano, columns, splits, labels)
            rows.append(
                dict(
                    folds=kind_of_folds,
                    map=NANO,
                    truth_from="",
                    predictors_from=h,
                    draw=-1,
                    **row,
                )
            )

        # the known maps: each from one half, analysed with the other's predictors
        mine = children[f * len(known) : (f + 1) * len(known)]
        for ((kind, truth_from), truth), child in zip(known.items(), mine):
            other = other_half(truth_from)
            columns = beyond_density.model_columns(cov[other], terms)
            rng = np.random.default_rng(child)
            for draw in range(n_noise):
                cohort = fake_cohort(truth, h_nano, rng)
                row = analyse(cohort, columns, splits, labels)
                rows.append(
                    dict(
                        folds=kind_of_folds,
                        map=kind,
                        truth_from=truth_from,
                        predictors_from=other,
                        draw=draw,
                        **row,
                    )
                )
    out = pd.DataFrame(rows)
    out.insert(1, "n_structures", len(cal.structures))
    return out


def left_summary(table: pd.DataFrame, kind: str = FLOOR_MAP, folds: str = RANDOM) -> dict:
    """The median and range of the share left, and of the replication, for one map.

    `table` is calibration.csv; with `kind` NANO, the rows of the two halves'
    predictors (not the merged one).
    """
    mine = table[(table["map"] == kind) & (table["folds"] == folds)]
    if kind == NANO:
        mine = mine[mine["predictors_from"].isin(HALVES)]
    return dict(
        left_median=float(mine["left"].median()),
        left_lo=float(mine["left"].min()),
        left_hi=float(mine["left"].max()),
        replication_median=float(mine["replication"].median()),
        replication_lo=float(mine["replication"].min()),
        replication_hi=float(mine["replication"].max()),
    )


def measured_once(split: list[str], terms: dict[str, tuple[str, ...]]) -> list[str]:
    """The main model's genes and measured densities that are the same in both halves.

    A gene of the model (its abundance genes, and the markers when the panel is in
    it) is the same when it has one usable experiment (not in `split`); a measured
    density always is, being one mouse measured once. Their mismatch is not in the
    floor, which errs low by it.
    """
    genes = list(terms["abundance"])
    if "markers" in terms["density"]:
        genes += list(beyond_density.MARKERS)
    once = [g for g in genes if g not in split]
    return once + [m for m in terms["density"] if m in synaptome.MEASURES]


# ===== The jackknife over structures =====


def left_out(n: int) -> int:
    """How many of `n` structures a subsample leaves out (jackknife_share of them)."""
    return int(round(BEYOND_CALIBRATION["jackknife_share"] * n))


def jackknife_subsamples(
    inputs: beyond_density.Inputs,
    splits: list[tuple[list[int], list[int]]],
    rng: np.random.Generator,
) -> Iterator[tuple[np.ndarray, beyond_density.Inputs, list]]:
    """The subsamples of the delete-d jackknife of the inputs' structures.

    beyond_calibration.n_jackknife of them, each leaving out left_out of the
    structures, drawn at random, then beyond_calibration.jackknife_splits of
    `splits`, drawn at random. Yields the positions kept, the inputs on them and
    the splits.
    """
    n = len(inputs.structures)
    d = left_out(n)
    for _ in range(BEYOND_CALIBRATION["n_jackknife"]):
        keep = np.sort(rng.choice(n, n - d, replace=False))
        sub = beyond_density.restrict(inputs, [inputs.structures[i] for i in keep])
        n_splits = BEYOND_CALIBRATION["jackknife_splits"]
        drawn = rng.choice(len(splits), n_splits, replace=False)
        yield keep, sub, [splits[i] for i in drawn]


def jackknife_table(rows: list[dict], n: int) -> pd.DataFrame:
    """The subsamples' rows as one table, the structures and those left out first."""
    out = pd.DataFrame(rows)
    out.insert(0, "n_left_out", left_out(n))
    out.insert(0, "n_structures", n)
    return out


def jackknife_sd(values: np.ndarray, n: int, d: int) -> float:
    """The delete-d jackknife SD of a statistic from its values on the subsamples.

    Each subsample leaves out d of n structures; the variance of the full-sample
    statistic is (n - d) / (d N) times the sum of squares about the subsamples'
    mean, N the number of subsamples.
    """
    values = np.asarray(values, float)
    centred = values - values.mean()
    return float(np.sqrt((n - d) / (d * len(values)) * np.sum(centred**2)))


def around(point: float, sd: float) -> tuple[float, float]:
    """The 95% interval of a value from its jackknife SD."""
    return float(point - 1.96 * sd), float(point + 1.96 * sd)


def jackknife_interval(
    table: pd.DataFrame, values: pd.Series, point: float
) -> tuple[float, float]:
    """The 95% interval of `point` from its `values` on the subsamples of `table`."""
    n, d = int(table["n_structures"].iloc[0]), int(table["n_left_out"].iloc[0])
    return around(point, jackknife_sd(values, n, d))


def paired_jackknife(
    inputs: beyond_density.Inputs,
    halves: dict[str, dict[str, dict[str, float]]],
    seed: int = 0,
) -> pd.DataFrame:
    """The nano map, the floor and the Gria1 map on the same subsamples of structures.

    On each subsample of the calibration structures (jackknife_subsamples): the
    predictors rebuilt from each half on the subsample, as production builds them;
    nano's share left by the main model with each half's predictors (their mean);
    the floor's and the Gria1 map's, each the mean over both directions and
    beyond_calibration.jackknife_draws draws of animals. Columns: nano, floor,
    gria1, nano_minus_floor, nano_minus_gria1. The animals' noise is drawn once for
    every structure and cut to each subsample, so the subsamples differ only in
    their structures.
    """
    terms = inputs.terms
    cal = beyond_density.restrict(
        inputs, calibration_structures(inputs.structures, halves)
    )
    n = len(cal.structures)
    rng = np.random.default_rng(seed)
    n_draws = BEYOND_CALIBRATION["jackknife_draws"]
    noise = rng.standard_normal((2 * len(HALVES), n_draws, len(beyond_density.ADULTS), n))
    rows = []
    for keep, sub, splits in jackknife_subsamples(cal, beyond_density.half_splits(), rng):
        nano = sub.nano
        cov = half_covariates(sub, halves)
        y = beyond_density.full_map(nano)
        h_nano = float(np.mean(beyond_density.ceiling(nano, splits)[0]))
        nano_left = np.mean(
            [
                analyse(
                    nano,
                    beyond_density.model_columns(cov[h], terms),
                    splits,
                    replication=False,
                )["left"]
                for h in HALVES
            ]
        )
        left = {FLOOR_MAP: [], GRIA1_MAP: []}
        known = known_maps(y, cov, terms)
        for m, ((kind, truth_from), truth) in enumerate(known.items()):
            columns = beyond_density.model_columns(cov[other_half(truth_from)], terms)
            for draw in range(n_draws):
                cohort = fake_cohort(truth, h_nano, noise=noise[m, draw][:, keep])
                left[kind].append(
                    analyse(cohort, columns, splits, replication=False)["left"]
                )
        floor_left = float(np.mean(left[FLOOR_MAP]))
        gria1_left = float(np.mean(left[GRIA1_MAP]))
        rows.append(
            dict(
                nano=float(nano_left),
                floor=floor_left,
                gria1=gria1_left,
                nano_minus_floor=float(nano_left - floor_left),
                nano_minus_gria1=float(nano_left - gria1_left),
            )
        )
    return jackknife_table(rows, n)


def jackknife_budget(
    inputs: beyond_density.Inputs,
    splits: list[tuple[list[int], list[int]]],
    rng: np.random.Generator,
) -> pd.DataFrame:
    """The budget of the main model on the subsamples of its structures.

    Per subsample (jackknife_subsamples): the ceiling on it and the cumulative share
    of each step of the main model (abundance, + density, + autofluorescence), each
    over that subsample's ceiling, with every predictor rebuilt on the structures
    kept.
    """
    rows = []
    for _, sub, few in jackknife_subsamples(inputs, splits, rng):
        _, explainable = beyond_density.ceiling(sub.nano, few)
        y = beyond_density.full_map(sub.nano)
        cov, _, _ = beyond_density.covariates_for(sub)
        steps = beyond_density.budget(y, cov, inputs.terms, explainable)
        rows.append(
            dict(
                ceiling=explainable,
                abundance=steps[0],
                abundance_density=steps[1],
                model=steps[2],
                left=1 - steps[2],
            )
        )
    return jackknife_table(rows, len(inputs.structures))


def jackknife_density_rows(
    inputs: beyond_density.Inputs,
    splits: list[tuple[list[int], list[int]]],
    rng: np.random.Generator,
) -> pd.DataFrame:
    """The check rows of DENSITY_ROWS on the same subsamples of their structures.

    On the structures where PSD95 is measured, each row refitted on every subsample
    (jackknife_subsamples), every predictor rebuilt on the structures kept: the
    share left (left_<key>) and the density's share alone (alone_<key>), over that
    subsample's ceiling; then the differences of DENSITY_CONTRASTS. The rows differ
    only in their density and meet the same subsamples, so their differences are
    paired.
    """
    measured = beyond_density.on_measured(inputs, beyond_density.MEASURED)
    rows = []
    for _, sub, few in jackknife_subsamples(measured, splits, rng):
        _, explainable = beyond_density.ceiling(sub.nano, few)
        y = beyond_density.full_map(sub.nano)
        cov, _, _ = beyond_density.covariates_for(sub)
        row = {}
        for key, density in DENSITY_ROWS.items():
            terms = beyond_density.model_terms(density)
            row[f"left_{key}"] = 1 - beyond_density.budget(y, cov, terms, explainable)[-1]
            columns = beyond_density.model_columns(cov, terms, ("density",))
            row[f"alone_{key}"] = beyond_density.cv_r2(y, columns) / explainable
        rows.append(row)
    out = jackknife_table(rows, len(measured.structures))
    for name, first, second, kind in DENSITY_CONTRASTS:
        out[name] = out[f"{kind}_{first}"] - out[f"{kind}_{second}"]
    return out


def density_contrasts(variants: pd.DataFrame, jack: pd.DataFrame) -> dict:
    """DENSITY_CONTRASTS: each difference between check rows, with its 95% interval.

    The point is the difference of the rows of variants.csv; the interval comes from
    the same difference on each subsample of jackknife_density_rows (`jack`).
    """
    rows = variants.set_index("key")
    column = {"left": "left", "alone": "density_alone"}
    out = {}
    for name, first, second, kind in DENSITY_CONTRASTS:
        point = float(rows.loc[first, column[kind]] - rows.loc[second, column[kind]])
        paired = jack[f"{kind}_{first}"] - jack[f"{kind}_{second}"]
        out[name] = (point, jackknife_interval(jack, paired, point))
    return out


# ===== Reading back =====


def load_calibration() -> pd.DataFrame:
    """The table run_beyond_calibration wrote."""
    if not CALIBRATION.exists():
        raise FileNotFoundError(
            f"{CALIBRATION} not found: run run_beyond_calibration.py first"
        )
    return pd.read_csv(CALIBRATION, keep_default_na=False, na_values=[""])


def load_jackknife() -> pd.DataFrame:
    """The paired jackknife run_beyond_calibration wrote."""
    if not JACKKNIFE.exists():
        raise FileNotFoundError(f"{JACKKNIFE} not found: run run_beyond_calibration.py")
    return pd.read_csv(JACKKNIFE)


def main() -> None:
    """Split the experiments, run the known maps and the real one, write the table."""
    beyond_density.OUT.mkdir(parents=True, exist_ok=True)
    inputs = beyond_density.load_inputs()
    halves, split = half_profiles()
    once = measured_once(split, inputs.terms)
    print(
        f"{len(split)} genes have two halves of experiments; of the main model's "
        f"genes and measures, these are the same in both: {', '.join(once)}"
    )
    table = calibration_rows(inputs, halves)
    table.to_csv(CALIBRATION, index=False)
    print(f"{int(table['n_structures'].iloc[0])} structures -> {CALIBRATION}")
    for folds in (RANDOM, BLOCKS):
        print(f"  folds: {folds}")
        mine = table[table["folds"] == folds]
        for r in mine[mine["map"] == NANO].itertuples():
            print(
                f"    nano, predictors from {r.predictors_from:4s}  left {r.left:6.1%}, "
                f"leftover replicates {r.replication:.3f}"
            )
        for kind in (FLOOR_MAP, GRIA1_MAP):
            known = mine[mine["map"] == kind]
            print(
                f"    {kind:22s} left {known['left'].median():6.1%} "
                f"({known['left'].min():.1%} to {known['left'].max():.1%}), "
                f"leftover replicates {known['replication'].median():.3f} "
                f"({known['replication'].min():.3f} to "
                f"{known['replication'].max():.3f}), {len(known)} rows"
            )

    # the nano map against the floor and the Gria1 map on the same subsamples
    jack = paired_jackknife(inputs, halves)
    jack.to_csv(JACKKNIFE, index=False)
    n, d = int(jack["n_structures"].iloc[0]), int(jack["n_left_out"].iloc[0])
    nano = left_summary(table, NANO)["left_median"]
    for other, label in ((FLOOR_MAP, "floor"), (GRIA1_MAP, "gria1")):
        point = nano - left_summary(table, other)["left_median"]
        lo, hi = jackknife_interval(jack, jack[f"nano_minus_{label}"], point)
        print(
            f"  nano minus {other}: {point:+.1%} (95% {lo:+.1%} to {hi:+.1%}, "
            f"jackknife of {len(jack)} subsamples leaving out {d} of {n})"
        )
