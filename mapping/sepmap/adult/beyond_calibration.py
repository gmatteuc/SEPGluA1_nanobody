"""The model of analysis 4 on maps whose answer is known: the calibration floor.

A leftover is what the model does not predict. Part of it can be real, and part is
the predictors measuring something slightly different from what the map is made of:
each Allen experiment is one P56 mouse on a 200 um grid, and two experiments of one
gene agree only in part (their reliability, ish.gene_table). So even a map made of
nothing but Gria1 mRNA and synapse density would leave something, because the Allen
maps that predict it are not the maps it is made of. The calibration measures how
much, with the main model of adult.beyond_density (Gria1 and synapse density,
straight), and with a check row's own model when adult.beyond_checks asks:

    two halves    each gene's Allen experiments (those the gene table uses) are
                  split in two, alternately by id (the 1st, 3rd, ... against the
                  2nd, 4th, ...), and each half is merged as the gene table merges
                  (ish.reliability.merge). Every gene of the main model has two or
                  more usable experiments, so each half holds its own; a gene
                  measured once would be the same in both halves and its mismatch
                  missing from the floor (measured_once lists any). The predictors
                  are built from each half (build_covariates), psd_pc1 from one list
                  of genes for both, those with two halves (shared_psd_genes);
                  autofluorescence, measured in our own brains, is the same in both,
                  and so is a measured density: one mouse, measured once
    known map     the model fitted to the nano map with one half's predictors, its
                  fitted values: exactly what Gria1 and synapse density predict,
                  with nano's own weights, and nothing else
    animals       ten made-up adults per map: the known map (ranked, scaled to unit
                  variance) plus noise independent per structure and adult, of
                  variance 5 (1 / h - 1), so that two halves of five agree as nano's
                  halves do (h, their mean agreement)
    analysis      the same model with the other half's predictors: the ceiling, the
                  CV R2, the share left, and how well the two half-cohort leftovers
                  agree; both ways round, with beyond_calibration.n_noise draws of
                  the animals each way. And with the same half's predictors (the
                  rows of OWN_HALF), which leaves only the animals' noise
    nano          the real map with each half's predictors, read on the same footing
                  as the known map, and with the merged predictors of the
                  production run (predictors_from "both") on the same structures
    structures    those of the model where both halves measure every gene of it, as
                  the fit's own rule asks of the merged profiles: a structure that
                  only one half measures drops out
    folds         the main model twice: with the production folds (random, averaged
                  over beyond.cv_repeats shufflings) and with folds of spatial blocks
                  (beyond_density.block_labels), so the fit with blocked folds has its
                  own floor

The known map gives the floor: the share any map would leave from one Allen map
disagreeing with another. Its leftover replicates across the made-up animals too,
which is why replication alone cannot tell biology from ISH mismatch.

The noise is independent between structures, where real animals deviate smoothly;
the floor calibrates a size, it does not model the cohort. Nor is the floor exact,
and its errors run both ways. What it misses makes it err low: the mismatch of
Allen's P56 mice with these brains (age, strain, the 200 um grid against 20 um
masks, registration), which nano meets and the known map does not. What it holds
makes it err high: the known map is built from one half of the experiments and read
with the other, so it meets the errors of two halves, where nano, read with one
half, meets one. Read with its own half the known map leaves only the animals'
noise, so the floor less that is the two halves' disagreement, of which one half's
error is about half; the floor errs low only where the P56 mismatch costs more.

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
it. The nano map and the floor meet the same subsamples of the calibration's
structures (paired_jackknife), so their difference has an interval of its own; the
partition and the weights get theirs on the subsamples of the main model's
structures (jackknife_partition).

Neighbouring structures are not independent: the map and its predictors vary
smoothly, and the folds of spatial blocks leave more than random folds. So every
interval is also given by a jackknife over spatial blocks (block_subsamples): each
subsample leaves out one block of beyond.cv_blocks, the k-means clusters of the
structures' centroids that the folds of spatial blocks use, and the SD is that of
the grouped jackknife (block_jackknife_sd). It is wider, and with twenty blocks
itself uncertain; both are reported, the random one as the method fixed first.

Writes, in adult_v2/ish_analysis/beyond/ under the data root:

    calibration.csv            one row per folds, map, direction and draw (draw -1
                               for the real map): the structures, the half
                               agreement, the ceiling, the CV R2, the share left and
                               the leftover's replication
    calibration_jackknife.csv  per subsample of the calibration's structures: the
                               nano map's share left, the floor's, their difference
    calibration_jackknife_blocks.csv  the same with one spatial block left out at a
                               time

Run by run_beyond_calibration.py; adult.beyond_checks runs the same on each check
row, adult.beyond_figures the jackknife of the partition.
"""

from collections.abc import Iterator

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from sepmap.adult import beyond_density, synaptome
from sepmap.config import SETTINGS
from sepmap.ish import gene_table, reliability
from sepmap.structures import centroid_xyz

# the draws of animal noise per direction; the jackknife's size
BEYOND_CALIBRATION = SETTINGS["beyond_calibration"]

CALIBRATION = beyond_density.OUT / "calibration.csv"
JACKKNIFE = beyond_density.OUT / "calibration_jackknife.csv"
JACKKNIFE_BLOCKS = beyond_density.OUT / "calibration_jackknife_blocks.csv"

# the maps, as calibration.csv names them: the real one, the floor (a map made of
# Gria1 and synapse density, read with the other half's predictors) and the same map
# read with its own half's, which leaves only the animals' noise
NANO = "nano"
FLOOR_MAP = "Gria1 and density"
OWN_HALF = "Gria1 and density, own half"

# the halves, as calibration.csv names them, and the predictors of the production
# run, both halves merged
HALVES = ("A", "B")
MERGED = "both"

# the two kinds of folds, as calibration.csv names them
RANDOM = "random"
BLOCKS = "spatial blocks"


# ===== The halves and the known map =====


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
    structures: list[str],
    halves: dict[str, dict[str, dict[str, float]]],
    genes: tuple[str, ...],
) -> list[str]:
    """The structures of the fit where both halves measure every one of `genes`."""
    return [
        s for s in structures if all(s in halves[h][g] for h in HALVES for g in genes)
    ]


def calibration_inputs(
    inputs: beyond_density.Inputs,
    halves: dict[str, dict[str, dict[str, float]]],
    terms: dict[str, tuple[str, ...]],
) -> beyond_density.Inputs:
    """The inputs on the structures where both halves measure the model's genes."""
    genes = beyond_density.model_genes(terms)
    return beyond_density.restrict(
        inputs, calibration_structures(inputs.structures, halves, genes)
    )


def shared_psd_genes(inputs: beyond_density.Inputs, halves: dict) -> list[str]:
    """The genes of psd_pc1 in the calibration: one list for both halves.

    The postsynaptic-density genes measured on every structure of the inputs in both
    halves whose two halves differ, that is with two or more usable experiments: a
    gene measured once would be the same in both halves, and its mismatch missing
    from the floor; a list chosen per half would differ between the halves.
    """
    measured = [
        set(beyond_density.psd_genes_measured(halves[h], inputs.role, inputs.structures))
        for h in HALVES
    ]
    a, b = (halves[h] for h in HALVES)
    both = set.intersection(*measured)
    return sorted(g for g in both if a[g] != b[g])


def half_covariates(
    inputs: beyond_density.Inputs, halves: dict
) -> dict[str, dict[str, np.ndarray]]:
    """The predictors built from each half of the experiments, on the inputs' structures.

    Autofluorescence and the measured densities are the same in both; psd_pc1 is
    built from the same genes in both (shared_psd_genes).
    """
    psd = shared_psd_genes(inputs, halves)
    return {
        h: beyond_density.build_covariates(
            halves[h],
            inputs.role,
            inputs.auto,
            inputs.structures,
            synapses=inputs.synapses,
            psd_genes=psd,
        )[0]
        for h in HALVES
    }


def known_maps(
    y: np.ndarray,
    cov: dict[str, dict[str, np.ndarray]],
    terms: dict[str, tuple[str, ...]],
    bend: bool = False,
) -> dict[str, np.ndarray]:
    """The known map from each half, {half the truth comes from: fitted values}.

    The model `terms` (bent when `bend`) fitted to the nano map `y` with that half's
    predictors: what they predict of it, and nothing else.
    """
    out = {}
    for h in HALVES:
        columns = beyond_density.model_columns(cov[h], terms, bend=bend)
        out[h] = y - beyond_density.residual(y, columns)
    return out


def nano_rows(
    cal: beyond_density.Inputs,
    cov: dict[str, dict[str, np.ndarray]],
    terms: dict[str, tuple[str, ...]],
    bend: bool,
    splits: list[tuple[list[int], list[int]]],
    folds: tuple[str, list[np.ndarray] | None],
    replication: bool,
) -> list[dict]:
    """The real map with each half's predictors and with the merged ones, one folds."""
    kind_of_folds, labels = folds
    merged, _, _ = beyond_density.covariates_for(cal)
    rows = []
    for h, covariates in [(h, cov[h]) for h in HALVES] + [(MERGED, merged)]:
        columns = beyond_density.model_columns(covariates, terms, bend=bend)
        row = analyse(cal.nano, columns, splits, labels, replication)
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
    return rows


def known_rows(
    known: dict[str, np.ndarray],
    cov: dict[str, dict[str, np.ndarray]],
    model: tuple[dict[str, tuple[str, ...]], bool],
    reading: tuple[str, list, list[np.ndarray] | None],
    children: list[np.random.SeedSequence],
    noise: tuple[float, int, bool],
) -> list[dict]:
    """The known map from each half, analysed n_noise times, one folds.

    `model` is the terms and whether they bend; `reading` is the map's name in the
    table, the splits and the folds' labels; `noise` the half agreement the animals
    are drawn to, the draws per direction and whether to compute the replication.
    With the map FLOOR_MAP the known map is read with the other half's predictors,
    with OWN_HALF with its own. Each direction draws from its own `children` entry.
    """
    terms, bend = model
    kind, splits, labels = reading
    h_nano, n_noise, replication = noise
    rows = []
    for (truth_from, truth), child in zip(known.items(), children):
        half = other_half(truth_from) if kind == FLOOR_MAP else truth_from
        columns = beyond_density.model_columns(cov[half], terms, bend=bend)
        rng = np.random.default_rng(child)
        for draw in range(n_noise):
            cohort = fake_cohort(truth, h_nano, rng)
            row = analyse(cohort, columns, splits, labels, replication)
            rows.append(
                dict(map=kind, truth_from=truth_from, predictors_from=half, draw=draw)
                | row
            )
    return rows


def calibration_rows(
    inputs: beyond_density.Inputs,
    halves: dict[str, dict[str, dict[str, float]]],
    terms: dict[str, tuple[str, ...]] | None = None,
    bend: bool = False,
    blocks: bool = True,
    replication: bool = True,
    own_half: bool = False,
    n_noise: int | None = None,
    seed: int = 0,
) -> pd.DataFrame:
    """calibration.csv: the real map and the known map, analysed with each half.

    The model `terms` (the main model by default; bent when `bend`) throughout, on
    the inputs' structures where both halves measure its genes. The known map is
    built from one half and analysed with the other half's predictors, both ways
    round, n_noise draws each, with random folds and, when `blocks`, with folds of
    spatial blocks; with `own_half` also with its own half's predictors, random
    folds. Every (folds, direction) has a generator of its own, so a draw does not
    depend on which others ran. Without `replication` the leftovers' agreement is
    not computed (NaN), which a check row does not need.
    """
    if terms is None:
        terms = inputs.terms
    if n_noise is None:
        n_noise = BEYOND_CALIBRATION["n_noise"]
    cal = calibration_inputs(inputs, halves, terms)
    cov = half_covariates(cal, halves)
    splits = beyond_density.half_splits()
    y = beyond_density.full_map(cal.nano)
    h_nano = float(np.mean(beyond_density.ceiling(cal.nano, splits)[0]))
    folds = [(RANDOM, None)]
    if blocks:
        folds.append((BLOCKS, beyond_density.block_labels(centroid_xyz(cal.structures))))
    known = known_maps(y, cov, terms, bend)

    # the generators of the own-half reading come after those of the floor, so
    # asking for it leaves the floor's draws as they are
    n = len(known)
    children = np.random.SeedSequence(seed).spawn((len(folds) + 1) * n)
    noise = (h_nano, n_noise, replication)
    rows = []
    for f, (kind_of_folds, labels) in enumerate(folds):
        rows += nano_rows(
            cal, cov, terms, bend, splits, (kind_of_folds, labels), replication
        )
        found = known_rows(
            known,
            cov,
            (terms, bend),
            (FLOOR_MAP, splits, labels),
            children[f * n : (f + 1) * n],
            noise,
        )
        rows += [dict(folds=kind_of_folds) | r for r in found]
    if own_half:
        found = known_rows(
            known,
            cov,
            (terms, bend),
            (OWN_HALF, splits, None),
            children[len(folds) * n :],
            noise,
        )
        rows += [dict(folds=RANDOM) | r for r in found]
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


def nano_left(table: pd.DataFrame, folds: str = RANDOM) -> float:
    """Nano's share left on the calibration's structures: the mean over the two halves."""
    mine = table[(table["map"] == NANO) & (table["folds"] == folds)]
    return float(mine.loc[mine["predictors_from"].isin(HALVES), "left"].mean())


def own_half_left(table: pd.DataFrame) -> float:
    """The median share the known map leaves read with its own half: the animals' noise.

    `table` is calibration.csv with the rows of OWN_HALF; NaN without them.
    """
    mine = table.loc[table["map"] == OWN_HALF, "left"]
    return float(mine.median()) if len(mine) else float("nan")


def measured_once(split: list[str], terms: dict[str, tuple[str, ...]]) -> list[str]:
    """A model's genes and measured densities that are the same in both halves.

    A gene of the model (model_genes) is the same when it has one usable experiment
    (not in `split`); a measured density always is, being one mouse measured once.
    Their mismatch is not in the floor, which errs low by it. psd_pc1 has none: the
    calibration builds it from genes with two halves (shared_psd_genes).
    """
    once = [g for g in beyond_density.model_genes(terms) if g not in split]
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


def block_subsamples(
    inputs: beyond_density.Inputs,
    splits: list[tuple[list[int], list[int]]],
    rng: np.random.Generator,
) -> Iterator[tuple[np.ndarray, beyond_density.Inputs, list]]:
    """The subsamples of the jackknife over spatial blocks of the inputs' structures.

    One per block of beyond_density.spatial_blocks (beyond.cv_blocks k-means
    clusters of the centroids), leaving that block out, then
    beyond_calibration.jackknife_splits of `splits`, drawn at random. Yields as
    jackknife_subsamples does.
    """
    block = beyond_density.spatial_blocks(centroid_xyz(inputs.structures))
    for b in np.unique(block):
        keep = np.flatnonzero(block != b)
        sub = beyond_density.restrict(inputs, [inputs.structures[i] for i in keep])
        n_splits = BEYOND_CALIBRATION["jackknife_splits"]
        drawn = rng.choice(len(splits), n_splits, replace=False)
        yield keep, sub, [splits[i] for i in drawn]


def subsamples(
    inputs: beyond_density.Inputs,
    splits: list[tuple[list[int], list[int]]],
    rng: np.random.Generator,
    blocks: bool,
) -> Iterator[tuple[np.ndarray, beyond_density.Inputs, list]]:
    """block_subsamples when `blocks`, else jackknife_subsamples."""
    if blocks:
        return block_subsamples(inputs, splits, rng)
    return jackknife_subsamples(inputs, splits, rng)


def jackknife_table(rows: list[dict], n: int, blocks: bool = False) -> pd.DataFrame:
    """The subsamples' rows as one table, the structures and those left out first.

    Over spatial blocks, n_left_out is NaN and n_blocks holds the number of blocks,
    which tells jackknife_interval which jackknife made the table.
    """
    out = pd.DataFrame(rows)
    if blocks:
        out.insert(0, "n_blocks", len(rows))
        out.insert(0, "n_left_out", np.nan)
    else:
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


def block_jackknife_sd(values: np.ndarray) -> float:
    """The grouped jackknife SD of a statistic from its values, one block left out each.

    With G blocks the variance is (G - 1) / G times the sum of squares about the
    mean of the G values; blocks of unequal size make it approximate.
    """
    values = np.asarray(values, float)
    g = len(values)
    centred = values - values.mean()
    return float(np.sqrt((g - 1) / g * np.sum(centred**2)))


def around(point: float, sd: float) -> tuple[float, float]:
    """The 95% interval of a value from its jackknife SD."""
    return float(point - 1.96 * sd), float(point + 1.96 * sd)


def jackknife_interval(
    table: pd.DataFrame, values: pd.Series, point: float
) -> tuple[float, float]:
    """The 95% interval of `point` from its `values` on the subsamples of `table`.

    Over spatial blocks when the table has n_blocks, else the delete-d jackknife.
    """
    if "n_blocks" in table.columns:
        return around(point, block_jackknife_sd(values))
    n, d = int(table["n_structures"].iloc[0]), int(table["n_left_out"].iloc[0])
    return around(point, jackknife_sd(values, n, d))


def subsample_floor(
    sub: beyond_density.Inputs,
    keep: np.ndarray,
    splits: list[tuple[list[int], list[int]]],
    halves: dict,
    model: tuple[dict[str, tuple[str, ...]], bool],
    noise: np.ndarray,
) -> dict:
    """Nano's share left and the floor's on one subsample of the calibration's structures.

    `model` is the terms and whether they bend; `noise` holds the animals' noise
    for every structure (directions x draws x adults x structures), cut here to the
    structures kept, so the subsamples differ only in their structures.
    """
    terms, bend = model
    cov = half_covariates(sub, halves)
    y = beyond_density.full_map(sub.nano)
    h_nano = float(np.mean(beyond_density.ceiling(sub.nano, splits)[0]))
    nano = np.mean(
        [
            analyse(
                sub.nano,
                beyond_density.model_columns(cov[h], terms, bend=bend),
                splits,
                replication=False,
            )["left"]
            for h in HALVES
        ]
    )
    floor = []
    known = known_maps(y, cov, terms, bend)
    for m, (truth_from, truth) in enumerate(known.items()):
        columns = beyond_density.model_columns(
            cov[other_half(truth_from)], terms, bend=bend
        )
        for draw in range(noise.shape[1]):
            cohort = fake_cohort(truth, h_nano, noise=noise[m, draw][:, keep])
            floor.append(analyse(cohort, columns, splits, replication=False)["left"])
    floor_left = float(np.mean(floor))
    return dict(
        nano=float(nano), floor=floor_left, nano_minus_floor=float(nano - floor_left)
    )


def paired_jackknife(
    inputs: beyond_density.Inputs,
    halves: dict[str, dict[str, dict[str, float]]],
    terms: dict[str, tuple[str, ...]] | None = None,
    bend: bool = False,
    seed: int = 0,
    blocks: bool = False,
) -> pd.DataFrame:
    """The nano map and the floor on the same subsamples of the calibration's structures.

    The model `terms` (the main model by default; bent when `bend`). On each
    subsample (jackknife_subsamples, or block_subsamples when `blocks`): the
    predictors rebuilt from each half on the subsample, as production builds them;
    nano's share left with each half's predictors (their mean); the floor's, the
    mean over both directions and beyond_calibration.jackknife_draws draws of
    animals. Columns: nano, floor, nano_minus_floor. The animals' noise is drawn
    once for every structure and cut to each subsample, so the subsamples differ
    only in their structures.
    """
    if terms is None:
        terms = inputs.terms
    cal = calibration_inputs(inputs, halves, terms)
    n = len(cal.structures)
    rng = np.random.default_rng(seed)
    n_draws = BEYOND_CALIBRATION["jackknife_draws"]
    noise = rng.standard_normal((len(HALVES), n_draws, len(beyond_density.ADULTS), n))
    rows = []
    splits = beyond_density.half_splits()
    for keep, sub, few in subsamples(cal, splits, rng, blocks):
        rows.append(subsample_floor(sub, keep, few, halves, (terms, bend), noise))
    return jackknife_table(rows, n, blocks)


def jackknife_partition(
    inputs: beyond_density.Inputs,
    splits: list[tuple[list[int], list[int]]],
    rng: np.random.Generator,
    blocks: bool = False,
) -> pd.DataFrame:
    """The main model's partition and weights on the subsamples of its structures.

    Per subsample (jackknife_subsamples, or block_subsamples when `blocks`): the
    ceiling on it, the shares of Gria1 alone, density alone and both
    (beyond_density.model_shares), the four parts (beyond_density.partition) and the
    weights (weight_<term>), with every predictor rebuilt on the structures kept.
    """
    rows = []
    for _, sub, few in subsamples(inputs, splits, rng, blocks):
        _, explainable = beyond_density.ceiling(sub.nano, few)
        y = beyond_density.full_map(sub.nano)
        cov, _, _ = beyond_density.covariates_for(sub)
        shares = beyond_density.model_shares(y, cov, inputs.terms, explainable)
        row = dict(ceiling=explainable, **shares)
        row.update(beyond_density.partition(shares))
        for term, weight in beyond_density.weights(y, cov, inputs.terms).items():
            row[f"weight_{term}"] = weight
        rows.append(row)
    return jackknife_table(rows, len(inputs.structures), blocks)


# ===== Reading back =====


def load_calibration() -> pd.DataFrame:
    """The table run_beyond_calibration wrote."""
    if not CALIBRATION.exists():
        raise FileNotFoundError(
            f"{CALIBRATION} not found: run run_beyond_calibration.py first"
        )
    return pd.read_csv(CALIBRATION, keep_default_na=False, na_values=[""])


def load_jackknife(blocks: bool = False) -> pd.DataFrame:
    """The paired jackknife run_beyond_calibration wrote; over spatial blocks if `blocks`.

    One subsample per row, as paired_jackknife returns it.
    """
    path = JACKKNIFE_BLOCKS if blocks else JACKKNIFE
    if not path.exists():
        raise FileNotFoundError(f"{path} not found: run run_beyond_calibration.py")
    return pd.read_csv(path)


def main() -> None:
    """Split the experiments, run the known map and the real one, write the tables."""
    beyond_density.OUT.mkdir(parents=True, exist_ok=True)
    inputs = beyond_density.load_inputs()
    halves, split = half_profiles()
    once = measured_once(split, inputs.terms)
    print(
        f"{len(split)} genes have two halves of experiments; of the main model's "
        f"genes ({', '.join(beyond_density.model_genes(inputs.terms))}), the same in "
        f"both: {', '.join(once) or 'none'}"
    )
    table = calibration_rows(inputs, halves, own_half=True)
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
        known = mine[mine["map"] == FLOOR_MAP]
        print(
            f"    the floor ({FLOOR_MAP})  left {known['left'].median():6.1%} "
            f"({known['left'].min():.1%} to {known['left'].max():.1%}), "
            f"leftover replicates {known['replication'].median():.3f} "
            f"({known['replication'].min():.3f} to "
            f"{known['replication'].max():.3f}), {len(known)} rows"
        )
    print(
        f"  the same map read with its own half: left {own_half_left(table):6.1%}, the "
        "animals' noise alone"
    )

    # the nano map against the floor on the same subsamples, random and in blocks
    jack = paired_jackknife(inputs, halves)
    jack.to_csv(JACKKNIFE, index=False)
    n, d = int(jack["n_structures"].iloc[0]), int(jack["n_left_out"].iloc[0])
    point = nano_left(table) - left_summary(table)["left_median"]
    lo, hi = jackknife_interval(jack, jack["nano_minus_floor"], point)
    print(
        f"  nano minus the floor: {point:+.1%} (95% {lo:+.1%} to {hi:+.1%}, "
        f"jackknife of {len(jack)} subsamples leaving out {d} of {n})"
    )
    jack = paired_jackknife(inputs, halves, blocks=True)
    jack.to_csv(JACKKNIFE_BLOCKS, index=False)
    lo, hi = jackknife_interval(jack, jack["nano_minus_floor"], point)
    print(
        f"  over spatial blocks: 95% {lo:+.1%} to {hi:+.1%} ({len(jack)} blocks, one "
        "left out at a time)"
    )
