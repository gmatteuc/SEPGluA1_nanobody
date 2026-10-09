"""What Gria1 expression and synapse density leave of the adult map (analysis 4).

The question: is the adult nano map just how much GluA1 a structure makes, or just
how many synapses it has? Ranking genes against the map cannot tell, because a map
of pure synapse density also puts the postsynaptic genes on top (analysis 3). So the
question is turned round: predict the map from Gria1 expression and synapse density,
and ask how much is left, whether the leftover is more than one Allen map
disagreeing with another, and what it looks like.

    nano map rank  ~  Gria1 rank + synapse-density rank  ->  leftover

This is the main model of settings.toml [beyond], fixed before it was run
(docs/ISH_ANALYSIS.md, section 4.1); the check rows beside it ([beyond.checks]) are
run by adult.beyond_checks, the controls by adult.beyond_controls. The steps, each
printing its numbers:

    0  the structures   the declared set (sepmap.structures), less the structures
                        where Gria1 or a density gene has no ISH value
    1  the ceiling      how reproducible the map is: nothing can predict variance
                        that is not there, so every R2 is read against it
    2  the partition    what Gria1 alone, synapse density alone and the two together
                        predict on structures the fit has not seen, split into what
                        only Gria1 predicts, what only density does, what they share
                        and what is left; and the weight of each term
    3  the leftover     whether it replicates across mice, and the value the ceiling
                        and the fit alone imply it would replicate at
    4  where it lives   which structures and divisions carry it and how consistently
                        across adults, and which genes' maps follow it, against the
                        leftover's own spatial null with the model removed from
                        every surrogate (leftover_genes)

Everything is one number per structure, on ranks across the structures:

    y            the ten adults' mean zref per structure (the declared reference)
    abundance    Gria1 (beyond.abundance). The mice carry SEP-GluA1, and Gria1
                 alone encodes GluA1; Gria2 to Gria4 encode partner subunits the
                 nanobody does not see, whose availability sets GluA1's assembly
                 and trafficking, so they belong to what the leftover may hold,
                 not to abundance (the four subunits are the check row
                 four_subunits)
    density      synapse density: the mean rank of the postsynaptic genes that
                 adult.density_markers chose by their agreement with the measured
                 PSD95 punctum density, by a rule that read no nano value and left
                 out every gene that places or regulates AMPA receptors
                 (density_markers.chosen). mRNA sits in somata, and a postsynaptic
                 gene's is made by the neurons that receive the synapses, the side
                 the receptor is on
    straight     each term enters once, as itself. A known map made of the two
                 terms, the calibration's floor, then has no curvature the model
                 could miss; the check row curved bends both (x, x^2, x^3) to show how
                 much of the leftover is curvature. It is a check, not a bound: terms
                 bent further take a little more (control E's fifth powers)
    CV R2        each structure predicted from a fit that never saw it: the
                 structures are shuffled, cut into five folds, and each fold
                 predicted from the other four; the share of variance missed is
                 averaged over beyond.cv_repeats such shufflings (seeded), since
                 one shuffling alone moves the leftover by several points.
                 In-sample R2 always grows with terms, so the held-out one is
                 quoted. Folds of spatial blocks (beyond.cv_blocks clusters of
                 neighbouring structures) are a variant (adult.beyond_controls)

The ceiling: the ten adults are split into two fives every possible way (126
splits), the two half-maps are correlated, and Spearman-Brown turns their mean
agreement r into the reliability of the ten-adult map,

    ceiling = 2 r / (1 + r)

A reliability is already a share of variance (true over observed), so the ceiling is
not squared again. Each model's share is its CV R2 over the ceiling. With A the share
of Gria1 alone, D that of density alone and AD that of the two together, the
reproducible map splits into

    Gria1 only     AD - D
    density only   AD - A
    shared         A + D - AD
    left           1 - AD

which add to 1. Held out, unlike in-sample, a part can come out negative: a term
that adds nothing to the other still costs a little through overfitting (a part
"only" below zero), and two terms that each predict little alone can predict more
together (a shared part below zero). A part is reported as computed, signed and never
clipped, and the figures draw a negative part leftwards from zero.

The weights: the straight model fitted on every structure with the map and both
terms z-scored, so the two coefficients share one unit, the SD of the map per SD of
the term, and say how much each term carries once the other is in.

The leftover is quoted as a range (adult.beyond_calibration, a jackknife over the
structures, and one over spatial blocks of them) beside the calibration floor: what
the same model leaves of a map that is exactly Gria1 and synapse density, measured
with other Allen experiments. Part of any leftover is one Allen map disagreeing with
another, and the floor says how much; the difference between the two is taken on the
same structures and resampled with them.

The leftover replicating across mice is not separate evidence. If the map
replicates, what is left of it once a smooth fit is removed must replicate too;
implied_replication gives the value that the ceiling and the fit alone predict, and
the figures show it beside the observed one. Nothing here measures what the
leftover is: it is what Gria1 expression and synapse density do not predict
(docs/ISH_ANALYSIS.md, sections 2 and 4.5, on the words and the reading).

The choices, and why:

    structures, not voxels  the gene data exist only at 200 um and in another mouse,
                            so at the voxel level the comparison would be mostly
                            registration error
    ranks everywhere        Allen expression energy has an arbitrary scale per
                            experiment; ranks assume only that more mRNA gives more
                            signal
    the declared set        grey matter seen in all ten adults: fibre tracts have no
                            synapses, and structures seen in a few mice do not
                            replicate
    naive and RWS pooled    the groups differ in whisker experience, not in what the
                            stain is; control D checks that the leftover does not
                            depend on the group

Writes, in adult_v2/ish_analysis/beyond/ under the data root:

    structures_used.csv        every structure of the adult table, in the fit or
                               not, why, and whether PSD95 is measured there
    variance_partition.csv     what each model predicts, against the ceiling
    partition.csv              the four parts of the reproducible map
    weights.csv                the main model's two weights
    replication.csv            per split of the adults, how well the two half-maps
                               and their two leftovers agree
    residual_by_structure.csv  per structure: map, prediction, leftover, and the
                               share of half-cohort leftovers with the same sign
    residual_by_division.csv   per division: the mean leftover, with its 95% over
                               resampled adults
    leftover_genes.csv        every gene's rho with the leftover, its spatial p
    leftover_sets.csv          the gene sets of analysis 3 against the leftover
    leftover_null.npz          the leftover's surrogates and every gene's rho with them

Run by run_beyond_density.py, which draws the working figures fig0 to fig3.
"""

import dataclasses
import itertools
from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.cluster.vq import kmeans2
from scipy.stats import rankdata, spearmanr

from sepmap import structures
from sepmap.adult import profiles, synaptome
from sepmap.adult.profiles import ADULTS
from sepmap.config import SETTINGS
from sepmap.ish import gene_ranking, gene_sets, gene_table, spatial_null
from sepmap.structures import load_centroids, load_structure_set

# the half-cohort size, the terms of the main model and of the check rows, the folds;
# the agreement a leftover needs to count as replicating; the density genes the rule
# of adult.density_markers chose
BEYOND = SETTINGS["beyond"]
DENSITY_MARKERS = SETTINGS["density_markers"]

OUT = structures.ISH_OUT / "beyond"
STRUCTURES_USED = OUT / "structures_used.csv"
PARTITION = OUT / "variance_partition.csv"
PARTS = OUT / "partition.csv"
WEIGHTS = OUT / "weights.csv"
REPLICATION = OUT / "replication.csv"
RESIDUALS = OUT / "residual_by_structure.csv"
RESIDUAL_DIVISIONS = OUT / "residual_by_division.csv"
LEFTOVER_GENES = OUT / "leftover_genes.csv"
LEFTOVER_SETS = OUT / "leftover_sets.csv"
LEFTOVER_NULL = OUT / "leftover_null.npz"

# the abundance term (Gria1); the four subunits of the check row four_subunits
ABUNDANCE = tuple(BEYOND["abundance"])
SUBUNITS = tuple(BEYOND["subunits"])

# the genes of each density composite, by the name of its predictor: the main model's
# (the genes the rule chose) and those of two check rows, the first proposal and the
# first version's marker panel
DENSITY_GENES = tuple(DENSITY_MARKERS["chosen"])
FIRST_PROPOSAL = tuple(BEYOND["first_proposal"])
MARKERS = tuple(BEYOND["markers"])
COMPOSITES = {
    "density": DENSITY_GENES,
    "first_proposal": FIRST_PROPOSAL,
    "markers": MARKERS,
}

# the groups of predictors, in the order a model lists them
ORDER = ("abundance", "density", "autofluorescence")

# the ontology panel's role whose genes make psd_pc1, the density component of a check
# row
PSD_ROLE = "control_psd"

# the parts of the reproducible map, in the order the tables and figures give them
PART_NAMES = ("gria1_only", "shared", "density_only", "left")


# ===== Loading =====


@dataclass
class Inputs:
    """What every step of analysis 4 reads, loaded once.

    `structures` are those the model runs on; `nano` and `auto` are adults x
    structures (zref of the nano and of the autofluorescence channel) in their
    order, and `synapses` the measured densities there (one column per
    synaptome.MEASURES, NaN where not measured). `expr` holds every gene's merged
    profile (mean ranks per structure), `role` its role in the ontology panel (''
    outside it); `structures_used` lists every structure of the adult table, in
    the fit or not, with the reason (structures_used.csv); `terms` is the main
    model, its predictors by group (model_terms).
    """

    structures: list[str]
    nano: np.ndarray
    auto: np.ndarray
    synapses: pd.DataFrame
    expr: dict[str, dict[str, float]]
    role: dict[str, str]
    p9_genes: set[str]
    division: dict[str, str]
    acronym: dict[str, str]
    structures_used: pd.DataFrame
    terms: dict[str, tuple[str, ...]]


def model_terms(
    density: Sequence[str] = ("density",),
    abundance: Sequence[str] = ABUNDANCE,
    autofluorescence: bool = False,
) -> dict[str, tuple[str, ...]]:
    """A model's predictors by group, named as build_covariates names them.

    Gria1 and the density composite by default, the main model; a check row changes
    the density or the abundance, or adds autofluorescence.
    """
    terms = {"abundance": tuple(abundance), "density": tuple(density)}
    if autofluorescence:
        terms["autofluorescence"] = ("autofluo",)
    return terms


def model_genes(terms: dict[str, tuple[str, ...]]) -> tuple[str, ...]:
    """The genes whose ISH values a model's predictors are made of.

    The abundance genes, and the genes of each density composite it uses; psd_pc1
    is made of whichever postsynaptic-density genes every structure has, and a
    measured density of no gene.
    """
    genes = list(terms["abundance"])
    for name in terms["density"]:
        genes += [g for g in COMPOSITES.get(name, ()) if g not in genes]
    return tuple(genes)


def structure_rows(
    set_table: pd.DataFrame,
    expr: dict[str, dict[str, float]],
    genes: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Every structure of the adult table: used by the fit or not, and why.

    A declared structure is used when every one of `genes` (the main model's by
    default) has an ISH value there: ish.min_voxels voxels in at least one of its
    experiments, as the gene table's profiles hold it. Columns: structure, acronym,
    division, used, reason, missing_genes.
    """
    if genes is None:
        genes = model_genes(model_terms())
    rows = []
    for r in set_table.itertuples():
        missing = [g for g in genes if r.structure not in expr.get(g, {})]
        if not r.in_set:
            reason = f"not in the declared set: {r.reason}"
            missing = []
        elif missing:
            reason = f"not measured: {' '.join(missing)}"
        else:
            reason = ""
        rows.append(
            dict(
                structure=r.structure,
                acronym=r.acronym,
                division=r.division,
                used=reason == "",
                reason=reason,
                missing_genes=" ".join(missing),
            )
        )
    return pd.DataFrame(rows)


def measured_structures(synapses: pd.DataFrame, measures: Sequence[str]) -> list[str]:
    """The structures of `synapses` where every one of `measures` is measured.

    `synapses` is indexed by structure.
    """
    values = synapses.reindex(columns=list(measures))
    return list(values.index[values.notna().all(axis=1)])


def load_inputs(genes: Sequence[str] | None = None, measured: bool = True) -> Inputs:
    """What every step reads: profiles, gene profiles, synapse density, structures.

    The structures are the declared ones where every one of `genes` is measured
    (structure_rows; the main model's genes by default, a check row's otherwise).
    With `measured` False the synaptome is not read: run_synaptome, which writes
    the synaptome's table, reads the fit this way.
    """
    per_mouse = profiles.load_per_mouse()
    set_table = load_structure_set()
    expr = gene_table.load_profiles()
    table = gene_table.per_gene(gene_table.load_gene_table())
    role = dict(zip(table["symbol"], table["ontology_role"]))
    p9_genes = set(table.loc[table["p9_gene"], "symbol"])
    used = structure_rows(set_table, expr, genes)
    fit = sorted(used.loc[used["used"], "structure"])

    # the measured densities on the fit, and where PSD95 is measured
    if measured:
        synapses = synaptome.load_density().reindex(fit)[list(synaptome.MEASURES)]
    else:
        synapses = pd.DataFrame(index=fit)
    covered = measured_structures(synapses, (synaptome.MEASURE,))
    used["psd95_measured"] = used["structure"].isin(covered)
    return Inputs(
        structures=fit,
        nano=gene_ranking.adult_matrix(per_mouse, "zref_nano", fit, ADULTS),
        auto=gene_ranking.adult_matrix(per_mouse, "zref_auto", fit, ADULTS),
        synapses=synapses,
        expr=expr,
        role=role,
        p9_genes=p9_genes,
        division=dict(zip(set_table["structure"], set_table["division"])),
        acronym=dict(zip(set_table["structure"], set_table["acronym"])),
        structures_used=used,
        terms=model_terms(),
    )


def restrict(inputs: Inputs, subset: Sequence[str]) -> Inputs:
    """The same inputs on some of their structures, in the order given."""
    columns = [inputs.structures.index(s) for s in subset]
    return dataclasses.replace(
        inputs,
        structures=list(subset),
        nano=inputs.nano[:, columns],
        auto=inputs.auto[:, columns],
        synapses=inputs.synapses.reindex(subset),
    )


def on_measured(inputs: Inputs, measures: Sequence[str]) -> Inputs:
    """The inputs on the structures where every one of `measures` is measured."""
    return restrict(inputs, measured_structures(inputs.synapses, measures))


# ===== Statistics =====


def composite(
    genes: Sequence[str], expr: dict[str, dict[str, float]], structures: list[str]
) -> np.ndarray:
    """One predictor from a set of genes: the rank of the mean of their rank profiles.

    Ranks rather than values because each Allen experiment carries its own
    arbitrary intensity scale, so raw numbers are not comparable between genes
    even when their orderings are; ranked again, so the composite is a rank across
    the structures as every other term.
    """
    mean = np.mean([rankdata([expr[g][s] for s in structures]) for g in genes], axis=0)
    return rankdata(mean)


def gene_matrix(
    expr: dict[str, dict[str, float]], structures: list[str], genes: Sequence[str]
) -> np.ndarray:
    """Rank profiles of many genes as one array, genes x structures.

    Ranks, since each Allen experiment has its own arbitrary intensity scale.
    """
    return np.array([rankdata([expr[g][s] for s in structures]) for g in genes])


def first_pc(
    genes: Sequence[str], expr: dict[str, dict[str, float]], structures: list[str]
) -> tuple[np.ndarray, float]:
    """The dominant shared axis of a gene set, and the share of variance it carries.

    Used for the postsynaptic-density genes: rather than naming a handful of
    markers and arguing about the choice, take whatever those genes have most in
    common and call that the density axis.
    """
    m = gene_matrix(expr, structures, genes)
    m = (m - m.mean(axis=1, keepdims=True)) / m.std(axis=1, keepdims=True)
    _, sv, vt = np.linalg.svd(m - m.mean(axis=0), full_matrices=False)
    pc = vt[0]

    # point it the same way as the genes
    if np.corrcoef(pc, m.mean(axis=0))[0, 1] < 0:
        pc = -pc
    return pc, float(sv[0] ** 2 / (sv**2).sum())


def projection(design: np.ndarray) -> np.ndarray:
    """The hat matrix of `design` (structures x columns).

    A map times its transpose is the map's least-squares fit on the columns.
    """
    return design @ np.linalg.pinv(design)


def residual(y: np.ndarray, predictors: Sequence[np.ndarray]) -> np.ndarray:
    """What is left of y after least squares on the predictors, plus an intercept."""
    design = np.column_stack(list(predictors) + [np.ones(len(y))])
    return y - design @ np.linalg.lstsq(design, y, rcond=None)[0]


def r_squared(y: np.ndarray, predictors: Sequence[np.ndarray]) -> float:
    """Variance explained on the data the fit was made from, so optimistic."""
    return float(1 - residual(y, predictors).var() / y.var())


def flexible(predictors: Sequence[np.ndarray]) -> list[np.ndarray]:
    """The same predictors, allowed to bend: x, x^2 and x^3 of each.

    The check row curved: a straight line through two rank variables assumes the
    relationship is not only monotone but evenly paced, and bending each term shows
    how much of the leftover is curvature the straight model misses.
    """
    return list(predictors) + [x**2 for x in predictors] + [x**3 for x in predictors]


def fold_labels(
    n: int, folds: int = 5, repeats: int | None = None, seed: int = 0
) -> list[np.ndarray]:
    """The fold of each of n structures, one array per shuffling.

    Each shuffling orders the structures at random and deals them into the folds in
    turn; the first, with seed 0, is the one shuffling of folds.csv.
    beyond.cv_repeats shufflings by default.
    """
    if repeats is None:
        repeats = BEYOND["cv_repeats"]
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(repeats):
        order = rng.permutation(n)
        label = np.empty(n, dtype=int)
        label[order] = np.arange(n) % folds
        out.append(label)
    return out


def spatial_blocks(
    xyz: np.ndarray, n_blocks: int | None = None, seed: int = 0
) -> np.ndarray:
    """The spatial block of each structure: k-means clusters of the centroids.

    `xyz` is structures x 3 (mm); beyond.cv_blocks clusters by default, seeded. The
    folds of spatial blocks and the jackknife over blocks both use them.
    """
    if n_blocks is None:
        n_blocks = BEYOND["cv_blocks"]
    _, block = kmeans2(np.asarray(xyz, float), n_blocks, seed=seed, minit="++")
    return block


def block_labels(
    xyz: np.ndarray,
    n_blocks: int | None = None,
    folds: int = 5,
    repeats: int | None = None,
    seed: int = 0,
) -> list[np.ndarray]:
    """Folds of spatial blocks: neighbouring structures held out together.

    The structures' centroids (xyz, structures x 3, mm) are cut into
    beyond.cv_blocks clusters (spatial_blocks); each shuffling deals the clusters
    into the folds in turn. A structure is then predicted from a fit that saw none
    of its neighbours, so smooth gradients cannot carry the prediction across.
    """
    if n_blocks is None:
        n_blocks = BEYOND["cv_blocks"]
    if repeats is None:
        repeats = BEYOND["cv_repeats"]
    block = spatial_blocks(xyz, n_blocks, seed)
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(repeats):
        fold_of_block = np.empty(n_blocks, dtype=int)
        fold_of_block[rng.permutation(n_blocks)] = np.arange(n_blocks) % folds
        out.append(fold_of_block[block])
    return out


def cv_predictions(
    y: np.ndarray,
    predictors: Sequence[np.ndarray],
    labels: list[np.ndarray] | None = None,
) -> np.ndarray:
    """Each structure's prediction from fits that never saw it, shufflings x structures.

    `labels` holds the fold of each structure per shuffling (fold_labels by default);
    each fold is predicted from a least-squares fit on the others.
    """
    n = len(y)
    if labels is None:
        labels = fold_labels(n)
    design = np.column_stack(list(predictors) + [np.ones(n)])
    out = np.empty((len(labels), n))
    for r, label in enumerate(labels):
        for k in np.unique(label):
            test = label == k
            beta = np.linalg.lstsq(design[~test], y[~test], rcond=None)[0]
            out[r, test] = design[test] @ beta
    return out


def cv_predict(
    y: np.ndarray,
    predictors: Sequence[np.ndarray],
    labels: list[np.ndarray] | None = None,
) -> np.ndarray:
    """Each structure's held-out prediction, the mean over the shufflings."""
    return cv_predictions(y, predictors, labels).mean(axis=0)


def cv_r2(
    y: np.ndarray,
    predictors: Sequence[np.ndarray],
    labels: list[np.ndarray] | None = None,
) -> float:
    """Variance explained on structures the fit has never seen (cv_predictions).

    The number to quote once a model has more than one term: adding predictors
    always improves the in-sample fit, and only held-out structures can say whether
    it improved the prediction. 1 minus the variance missed, averaged over the
    shufflings, over the variance of y.
    """
    predicted = cv_predictions(y, predictors, labels)
    missed = np.var(y[None, :] - predicted, axis=1).mean()
    return float(1 - missed / np.var(y))


def half_splits() -> list[tuple[list[int], list[int]]]:
    """Every way of cutting ten animals into two fives, each split counted once.

    A split is counted once by keeping only the halves that hold animal 0.
    """
    return [
        (list(p), [i for i in range(len(ADULTS)) if i not in p])
        for p in itertools.combinations(range(len(ADULTS)), BEYOND["half"])
        if 0 in p
    ]


def spearman_brown(r: float) -> float:
    """Reliability of a whole cohort from the agreement of its two halves; NaN at -1."""
    return 2 * r / (1 + r) if r > -1 else float("nan")


def replicates(agreement: float) -> bool:
    """Whether two half-cohort maps or leftovers agree at beyond.replication."""
    return agreement >= BEYOND["replication"]


def half_map(matrix: np.ndarray, rows: Sequence[int]) -> np.ndarray:
    """The ranked mean map of the adults in `rows` (matrix: adults x structures)."""
    return rankdata(matrix[list(rows)].mean(axis=0))


def full_map(matrix: np.ndarray) -> np.ndarray:
    """The ranked mean map of every adult."""
    return half_map(matrix, range(matrix.shape[0]))


def ceiling(
    matrix: np.ndarray, splits: list[tuple[list[int], list[int]]]
) -> tuple[list[float], float]:
    """The two half-maps' agreement over every split, and the ceiling.

    The ceiling is the Spearman-Brown reliability of the mean agreement: the share
    of the ten-adult map's variance that is reproducible.
    """
    agreement = [
        float(spearmanr(half_map(matrix, a), half_map(matrix, b)).statistic)
        for a, b in splits
    ]
    return agreement, spearman_brown(float(np.mean(agreement)))


def leftover_agreement(
    matrix: np.ndarray,
    columns: Sequence[np.ndarray],
    splits: list[tuple[list[int], list[int]]],
) -> list[float]:
    """How well the leftover of one half-cohort matches the other's, per split.

    Each half's ranked map is fitted to the model's columns on its own, and the two
    leftovers are compared by Spearman.
    """
    return [
        float(
            spearmanr(
                residual(half_map(matrix, a), columns),
                residual(half_map(matrix, b), columns),
            ).statistic
        )
        for a, b in splits
    ]


def implied_replication(half_agreement: float, r2: float, terms: int, n: int) -> float:
    """The leftover's replication that the map's reliability and the fit alone imply.

    A half-map is signal (a share h, the half-against-half agreement) plus noise
    (1 - h). The fit removes a share r2 of the variance, all of it signal, and
    absorbs terms / n of the noise (terms counts the intercept). Two half leftovers
    then share the signal left, h - r2, out of h - r2 plus the noise left, so

        implied = (h - r2) / ((h - r2) + (1 - h) (1 - terms / n))

    A leftover that replicates near this value carries no evidence beyond the
    map's own reliability; only the ceiling and the calibration speak to it.
    """
    signal = half_agreement - r2
    noise = (1 - half_agreement) * (1 - terms / n)
    return float(signal / (signal + noise)) if signal > 0 else float("nan")


# ===== The model =====


def psd_genes_measured(
    expr: dict[str, dict[str, float]], role: dict[str, str], structures: list[str]
) -> list[str]:
    """The panel's postsynaptic-density genes measured on every one of `structures`."""
    return sorted(
        g
        for g in expr
        if role.get(g, "") == PSD_ROLE and all(s in expr[g] for s in structures)
    )


def build_covariates(
    expr: dict[str, dict[str, float]],
    role: dict[str, str],
    auto: np.ndarray,
    structures: list[str],
    synapses: pd.DataFrame | None = None,
    psd_genes: Sequence[str] | None = None,
) -> tuple[dict[str, np.ndarray], list[str], float]:
    """Every predictor by name, built in one place so every module builds them alike.

    Each a rank across `structures`: the four subunits, each density composite of
    COMPOSITES, psd_pc1, autofluo (the mean of `auto`, adults x structures) and each
    measured density of `synapses` (structures x measures). A predictor whose genes
    or measure are not there for every structure is left out, and a model that
    needs it stops (predictors). psd_pc1 is the first component of `psd_genes`, by
    default every postsynaptic-density gene measured on the structures
    (psd_genes_measured); the calibration gives both halves one list. Returns the
    predictors, the genes behind psd_pc1 and the share of their variance it carries.
    """
    out = {}
    for g in SUBUNITS:
        if all(s in expr.get(g, {}) for s in structures):
            out[g] = rankdata([expr[g][s] for s in structures])
    for name, genes in COMPOSITES.items():
        if all(s in expr.get(g, {}) for g in genes for s in structures):
            out[name] = composite(genes, expr, structures)
    if psd_genes is None:
        psd_genes = psd_genes_measured(expr, role, structures)
    psd = list(psd_genes)
    psd_pc1, share = first_pc(psd, expr, structures)
    out["psd_pc1"] = rankdata(psd_pc1)
    out["autofluo"] = rankdata(auto.mean(axis=0))
    if synapses is not None:
        for measure in synapses.columns:
            values = synapses[measure].reindex(structures).to_numpy(float)
            if np.isfinite(values).all():
                out[measure] = rankdata(values)
    return out, psd, share


def covariates_for(inputs: Inputs) -> tuple[dict[str, np.ndarray], list[str], float]:
    """build_covariates on the inputs' own structures, measured densities included."""
    return build_covariates(
        inputs.expr, inputs.role, inputs.auto, inputs.structures, inputs.synapses
    )


def predictors(
    covariates: dict[str, np.ndarray],
    terms: dict[str, tuple[str, ...]],
    groups: Sequence[str] | None = None,
) -> list[np.ndarray]:
    """The straight predictors of the groups named (every group of `terms` by default).

    A predictor that build_covariates could not make on these structures stops the
    run: the structures were chosen for another model.
    """
    if groups is None:
        groups = [g for g in ORDER if g in terms]
    out = []
    for group in groups:
        for name in terms[group]:
            if name not in covariates:
                raise ValueError(
                    f"the predictor {name} is not measured on every structure of these "
                    "inputs: load them with that model's genes (model_genes)"
                )
            out.append(covariates[name])
    return out


def model_columns(
    covariates: dict[str, np.ndarray],
    terms: dict[str, tuple[str, ...]],
    groups: Sequence[str] | None = None,
    bend: bool = False,
) -> list[np.ndarray]:
    """The columns of a model: the groups' predictors, straight unless `bend`.

    model_columns(covariates, inputs.terms) is the main model: Gria1 and synapse
    density, each once.
    """
    straight = predictors(covariates, terms, groups)
    return flexible(straight) if bend else straight


def model_shares(
    y: np.ndarray,
    covariates: dict[str, np.ndarray],
    terms: dict[str, tuple[str, ...]],
    explainable: float,
    labels: list[np.ndarray] | None = None,
    bend: bool = False,
) -> dict[str, float]:
    """A model's shares of the reproducible map, held out: abundance, density, both.

    The CV R2 over the ceiling of the abundance terms alone, of the density terms
    alone and of the whole model (with autofluorescence when the model has it), and
    the share the whole model leaves (folds `labels`, fold_labels by default).
    """
    out = {}
    for key, groups in (
        ("abundance", ("abundance",)),
        ("density", ("density",)),
        ("model", None),
    ):
        columns = model_columns(covariates, terms, groups, bend)
        out[key] = cv_r2(y, columns, labels) / explainable
    out["left"] = 1 - out["model"]
    return out


def partition(shares: dict[str, float]) -> dict[str, float]:
    """The four parts of the reproducible map from the shares of a two-group model.

    `shares` holds abundance (A), density (D) and model (AD), as model_shares gives
    them: Gria1 only AD - D, density only AD - A, shared A + D - AD, left 1 - AD.
    They add to 1; a part may come out negative held out, and is kept as it is.
    """
    a, d, ad = shares["abundance"], shares["density"], shares["model"]
    return dict(gria1_only=ad - d, shared=a + d - ad, density_only=ad - a, left=1 - ad)


def weights(
    y: np.ndarray, covariates: dict[str, np.ndarray], terms: dict[str, tuple[str, ...]]
) -> dict[str, float]:
    """The straight model's coefficients with the map and every term z-scored.

    Fitted on every structure. Each is the SD of the map one SD of the term moves
    once the other terms are in, so the terms' weights compare directly.
    """
    names = [name for group in ORDER if group in terms for name in terms[group]]
    columns = predictors(covariates, terms)
    z = [(x - x.mean()) / x.std() for x in columns]
    target = (y - y.mean()) / y.std()
    beta = np.linalg.lstsq(np.column_stack(z), target, rcond=None)[0]
    return dict(zip(names, (float(b) for b in beta)))


# ===== Reading back =====


def load_leftover_genes() -> pd.DataFrame:
    """leftover_genes.csv as step 4 of this module wrote it, empty text where blank."""
    if not LEFTOVER_GENES.exists():
        raise FileNotFoundError(
            f"{LEFTOVER_GENES} not found: run run_beyond_density.py first"
        )
    genes = pd.read_csv(LEFTOVER_GENES, keep_default_na=False, na_values=[""])
    for column in ("in_model", "gene_sets"):
        genes[column] = genes[column].fillna("")
    return genes


def load_leftover_null(structures: list[str]) -> dict:
    """leftover_null.npz: the leftover's surrogates and every gene's rho with them.

    Returns structures, surrogates (as floats), rho (genes x surrogates) and genes;
    the run stops unless they were drawn on `structures`, those of the main model.
    """
    with np.load(LEFTOVER_NULL) as z:
        null = dict(
            structures=[str(s) for s in z["structures"]],
            surrogates=z["surrogates"].astype(float),
            rho=z["null_rho"],
            genes=[str(g) for g in z["genes"]],
        )
    if null["structures"] != structures:
        raise ValueError(
            f"{LEFTOVER_NULL} was drawn on other structures than the main model's "
            f"{len(structures)}; run run_beyond_density.py first"
        )
    return null


def load_partition() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """variance_partition.csv by key, partition.csv by part, weights.csv by term."""
    for path in (PARTITION, PARTS, WEIGHTS):
        if not path.exists():
            raise FileNotFoundError(f"{path} not found: run run_beyond_density.py first")
    return (
        pd.read_csv(PARTITION).set_index("key"),
        pd.read_csv(PARTS).set_index("part"),
        pd.read_csv(WEIGHTS).set_index("term"),
    )


# ===== Steps =====


def step0_structures(inputs: Inputs) -> None:
    """Write the structures the fit runs on, and why the others are left out."""
    print("\nstep 0  which structures the fit runs on")
    used = inputs.structures_used
    used.to_csv(STRUCTURES_USED, index=False)
    left_out = used[~used["used"]]
    n_fit = int(used["used"].sum())
    print(f"  {len(used)} structures in the adult table, {n_fit} in the fit")
    for reason, n in left_out["reason"].value_counts().items():
        print(f"    {n:3d}  {reason}")
    lacking = defaultdict(int)
    for text in left_out["missing_genes"]:
        for g in text.split():
            lacking[g] += 1
    if lacking:
        print(
            "  genes without a value in a declared structure: "
            + ", ".join(f"{g} {n}" for g, n in sorted(lacking.items()))
        )
    n_measured = int(used["psd95_measured"].sum())
    print(
        f"  synapse density: the mean rank of {', '.join(DENSITY_GENES)} "
        f"(density_markers.chosen); PSD95 puncta measured in {n_measured} of the "
        f"{n_fit}"
    )


def step1_ceiling(
    nano: np.ndarray, splits: list[tuple[list[int], list[int]]]
) -> tuple[list[float], float]:
    """How reproducible the map itself is: nothing below can beat this.

    Returns the two half-maps' agreement per split and the ceiling, as ceiling does.
    """
    print("\nstep 1  the ceiling: how much of this map is reproducible")
    agreement, explainable = ceiling(nano, splits)
    half = float(np.mean(agreement))
    print(f"  over {len(splits)} five-against-five splits of the ten adults:")
    print(f"    half-cohort against half-cohort   rho = {half:.3f}")
    print(f"    Spearman-Brown, all ten adults    {explainable:.3f}")
    print(f"  so {explainable:.1%} of the map's variance is reproducible, and every")
    print("  R2 below is read against that rather than against 1.")
    return agreement, explainable


def partition_models(
    covariates: dict[str, np.ndarray], terms: dict[str, tuple[str, ...]]
) -> list[tuple[str, str, list, bool]]:
    """The models of variance_partition.csv as (key, label, columns, bent).

    "abundance", "density" and "model" make the partition; the others describe it:
    autofluorescence alone, straight, and the main model curved (the check row
    curved). The four subunits are measured on fewer structures, so they are a check
    row of their own (adult.beyond_checks).
    """
    abundance = ", ".join(terms["abundance"])
    genes = ", ".join(DENSITY_GENES)
    auto = {"abundance": (), "density": (), "autofluorescence": ("autofluo",)}
    return [
        (
            "abundance",
            f"{abundance} alone",
            model_columns(covariates, terms, ("abundance",)),
            False,
        ),
        (
            "density",
            f"synapse density alone ({genes})",
            model_columns(covariates, terms, ("density",)),
            False,
        ),
        (
            "model",
            f"the main model: {abundance} + synapse density",
            model_columns(covariates, terms),
            False,
        ),
        (
            "autofluorescence",
            "autofluorescence alone",
            model_columns(covariates, auto, ("autofluorescence",)),
            False,
        ),
        (
            "curved",
            "the main model curved (x, x^2, x^3 of each term)",
            model_columns(covariates, terms, bend=True),
            True,
        ),
    ]


def partition_table(
    y: np.ndarray, models: list[tuple[str, str, list, bool]], explainable: float
) -> pd.DataFrame:
    """variance_partition.csv: each model's terms, R2 held out and in-sample, shares."""
    rows = []
    for key, label, columns, bent in models:
        cv = cv_r2(y, columns)
        rows.append(
            dict(
                key=key,
                model=label,
                terms=len(columns) + 1,
                bent=bent,
                cv_r2=cv,
                in_sample_r2=r_squared(y, columns),
                share_of_ceiling=cv / explainable,
                left=1 - cv / explainable,
            )
        )
    return pd.DataFrame(rows)


def parts_table(table: pd.DataFrame) -> pd.DataFrame:
    """partition.csv: the four parts of the reproducible map, from variance_partition."""
    shares = table.set_index("key")["share_of_ceiling"]
    parts = partition(shares)
    return pd.DataFrame(dict(part=list(parts), share=list(parts.values())))


def weights_table(
    y: np.ndarray, covariates: dict[str, np.ndarray], terms: dict[str, tuple[str, ...]]
) -> pd.DataFrame:
    """weights.csv: the main model's weights, and each term's rho with the map."""
    found = weights(y, covariates, terms)
    return pd.DataFrame(
        dict(
            term=list(found),
            weight=list(found.values()),
            rho_with_map=[float(spearmanr(y, covariates[t]).statistic) for t in found],
        )
    )


def step2_partition(
    y: np.ndarray,
    covariates: dict[str, np.ndarray],
    terms: dict[str, tuple[str, ...]],
    explainable: float,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """What Gria1 and synapse density predict, the partition and the weights.

    Returns variance_partition.csv's, partition.csv's and weights.csv's tables.
    """
    print("\nstep 2  what Gria1 and synapse density predict")
    print("  each predictor on its own (Spearman with the map):")
    for key in covariates:
        print(f"    {key:20s} rho {spearmanr(y, covariates[key]).statistic:+.3f}")

    table = partition_table(y, partition_models(covariates, terms), explainable)
    table.to_csv(PARTITION, index=False)
    print("  the models, scored on structures the fit has not seen:")
    for r in table.itertuples():
        print(
            f"    {r.model:52s} CV R2 {r.cv_r2:+.3f}   in-sample {r.in_sample_r2:.3f}"
            f"   {r.share_of_ceiling:6.1%} of the ceiling"
        )
    parts = parts_table(table)
    parts.to_csv(PARTS, index=False)
    print(
        "  the reproducible map: "
        + ", ".join(f"{r.part} {r.share:+.1%}" for r in parts.itertuples())
    )
    found = weights_table(y, covariates, terms)
    found.to_csv(WEIGHTS, index=False)
    print(
        "  the weights (z-scored): "
        + ", ".join(f"{r.term} {r.weight:+.3f}" for r in found.itertuples())
    )
    return table, parts, found


def replication_table(
    splits: list[tuple[list[int], list[int]]],
    map_agreement: list[float],
    left_agreement: list[float],
) -> pd.DataFrame:
    """replication.csv: per split, its two halves and the two agreements."""
    return pd.DataFrame(
        dict(
            half_a=[" ".join(ADULTS[i] for i in a) for a, _ in splits],
            half_b=[" ".join(ADULTS[i] for i in b) for _, b in splits],
            map_agreement=map_agreement,
            leftover_agreement=left_agreement,
        )
    )


def step3_replication(
    nano: np.ndarray,
    columns: list[np.ndarray],
    splits: list[tuple[list[int], list[int]]],
    map_agreement: list[float],
) -> tuple[list[float], float]:
    """Whether the leftover replicates across mice, and what that is worth."""
    print("\nstep 3  does the leftover replicate across mice?")
    agreement = leftover_agreement(nano, columns, splits)
    y = full_map(nano)
    h = float(np.mean(map_agreement))
    implied = implied_replication(h, r_squared(y, columns), len(columns) + 1, len(y))
    replication_table(splits, map_agreement, agreement).to_csv(REPLICATION, index=False)
    print(
        f"  half against half: the map {h:.3f}, its leftover {np.mean(agreement):.3f} "
        f"(lowest split {min(agreement):.3f})"
    )
    print(
        f"  the map's reliability and the fit alone imply {implied:.3f}: a leftover "
        "replicating near that carries no extra evidence"
    )
    return agreement, implied


def same_sign_share(
    nano: np.ndarray,
    columns: list[np.ndarray],
    splits: list[tuple[list[int], list[int]]],
    leftover: np.ndarray,
) -> np.ndarray:
    """Per structure, the share of the half-cohort leftovers with the cohort's sign.

    Every split gives two half-cohorts, each fitted on its own; a structure whose
    leftover has the same sign in all of them is above (or below) prediction in
    every five adults.
    """
    same = np.zeros(len(leftover))
    for a, b in splits:
        for half in (a, b):
            same += np.sign(residual(half_map(nano, half), columns)) == np.sign(leftover)
    return same / (2 * len(splits))


def residual_table(
    inputs: Inputs,
    y: np.ndarray,
    predicted: np.ndarray,
    leftover: np.ndarray,
    same: np.ndarray,
    per_adult: np.ndarray,
) -> pd.DataFrame:
    """residual_by_structure.csv: per structure, largest leftover first.

    Beside the cohort's leftover (residual): the share of half-cohort leftovers with
    its sign (same_sign), and the leftover of each adult's own map (`per_adult`,
    adults x structures) as a mean, an SD and a t across the adults (t_adults), the
    reliability that the bars of figure 04 show in grey.
    """
    n = per_adult.shape[0]
    mean = per_adult.mean(axis=0)
    sd = per_adult.std(axis=0, ddof=1)
    s = inputs.structures
    table = pd.DataFrame(
        dict(
            structure=s,
            acronym=[inputs.acronym.get(x, "") for x in s],
            division=[inputs.division.get(x, "") for x in s],
            nano_rank=y,
            predicted_rank=predicted,
            residual=leftover,
            same_sign=same,
            adult_mean=mean,
            adult_sd=sd,
            t_adults=mean / (sd / np.sqrt(n)),
        )
    )
    return table.sort_values("residual", ascending=False, ignore_index=True)


def per_adult_leftovers(nano: np.ndarray, columns: list[np.ndarray]) -> np.ndarray:
    """Each adult's own leftover, its ranked map fitted to the model.

    Adults x structures.
    """
    return np.array([residual(rankdata(row), columns) for row in nano])


def leftover_boot(
    nano: np.ndarray, columns: list[np.ndarray], seed: int = 0
) -> np.ndarray:
    """Leftovers of resampled cohorts, n_boot_mice x structures.

    Each is the mean of the adults drawn with replacement, ranked and fitted to the
    model's columns.
    """
    boot = gene_ranking.bootstrap_maps(nano, seed=seed)
    return np.array([residual(rankdata(b), columns) for b in boot])


def between_share(values: np.ndarray, groups: Sequence[str]) -> float:
    """The share of the variance of `values` that lies between `groups`.

    The variance of the group means, each weighted by its size, over the total: 0 when
    every group has the same mean, 1 when every value equals its group's mean.
    """
    values = np.asarray(values, float)
    frame = pd.DataFrame(dict(value=values, group=list(groups)))
    means = frame.groupby("group")["value"].transform("mean").to_numpy()
    return float(
        np.sum((means - values.mean()) ** 2) / np.sum((values - values.mean()) ** 2)
    )


def division_table(
    inputs: Inputs, leftover: np.ndarray, boot: np.ndarray
) -> pd.DataFrame:
    """residual_by_division.csv: per division, the mean leftover and its spread.

    The mean and median of the cohort's leftover over the division's structures, and
    the 2.5 and 97.5 percentiles of the mean over the cohorts of resampled adults
    (`boot`, resamples x structures, leftover_boot); highest mean first. A division
    whose interval excludes zero sits above or below prediction as a whole.
    """
    division = np.array([inputs.division.get(s, "") for s in inputs.structures])
    rows = []
    for name in sorted(set(division)):
        mine = division == name
        resampled = boot[:, mine].mean(axis=1)
        rows.append(
            dict(
                division=name,
                n_structures=int(mine.sum()),
                mean=float(leftover[mine].mean()),
                median=float(np.median(leftover[mine])),
                lo=float(np.percentile(resampled, 2.5)),
                hi=float(np.percentile(resampled, 97.5)),
            )
        )
    table = pd.DataFrame(rows)
    return table.sort_values("mean", ascending=False, ignore_index=True)


def model_membership(terms: dict[str, tuple[str, ...]]) -> dict[str, str]:
    """{gene: its part of the model}: 'abundance' for a term, 'density' in a composite."""
    out = {g: "abundance" for g in terms["abundance"]}
    for name in terms["density"]:
        out.update({g: "density" for g in COMPOSITES.get(name, ())})
    return out


def leftover_genes(
    inputs: Inputs,
    leftover: np.ndarray,
    columns: list[np.ndarray],
    boot: np.ndarray | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray, np.ndarray, list[str]]:
    """Every gene against the leftover, with the leftover's own spatial null.

    The leftover is a map of its own, with its own smoothness, so a gene's rho with
    it is tested as in analysis 1: against the leftover's surrogates (Burt 2020, on
    the structures' one-hemisphere centroids). The leftover is what a fit leaves,
    so it carries nothing of the model's columns; a surrogate drawn from it does,
    and against it a gene sharing the model's pattern would meet a null far too
    wide. So each surrogate goes through the same fit and only its residual is kept
    (a Freedman-Lane null with spatial surrogates). Gria1 enters the model directly
    and correlates with the leftover near zero by construction; the density genes
    enter only through their mean rank, and need not (in_model says which). A gene
    that follows the leftover predicts what Gria1 and synapse density leave. The
    gene sets of analysis 3 are read the same way. Cacng8 and the AMPA receptor
    complex family were named for the leftover in advance (ish.gene_sets); these
    tables describe every gene alike.

    `boot` holds the leftovers of resampled cohorts (leftover_boot, computed when
    not given). Returns the gene table, the set table, the surrogates, every gene's
    null rho (genes x surrogates) and the genes in the order of its rows.
    """
    s = inputs.structures
    centroids = load_centroids().loc[s]
    d = spatial_null.distance_matrix(centroids)
    surr = spatial_null.surrogates(leftover, d, seed=0)
    design = np.column_stack(list(columns) + [np.ones(len(leftover))])
    surr = surr - surr @ projection(design).T
    if boot is None:
        boot = leftover_boot(inputs.nano, columns)
    vectors = gene_ranking.gene_vectors(inputs.expr, s)
    table, null = gene_ranking.rank_genes(leftover, surr, boot, vectors)
    table = gene_ranking.with_q_and_ranks(table, inputs.p9_genes)
    table["in_model"] = table["symbol"].map(model_membership(inputs.terms)).fillna("")
    table["p9_gene"] = table["symbol"].isin(inputs.p9_genes)
    genes = gene_table.per_gene(gene_table.load_gene_table()).set_index("symbol")
    table["gene_sets"] = table["symbol"].map(genes["gene_sets"]).fillna("")
    table["reliability"] = table["symbol"].map(genes["reliability"])

    # the gene sets, as analysis 3 reads them on the map
    members = gene_sets.members_from_table(genes.reset_index())
    tested = list(vectors)
    rho = table.set_index("symbol")["rho"]
    sets = gene_sets.set_tests(members, {"leftover": rho}, {"leftover": null}, tested)
    table = table.sort_values("rho", ascending=False, ignore_index=True)
    return table, sets, surr, null, tested


def print_leftover_rows(table: pd.DataFrame) -> None:
    """Print rows of residual_by_structure.csv.

    Leftover in ranks, same-sign share, t across the adults' own leftovers.
    """
    for r in table.itertuples():
        print(
            f"    {r.structure[:48]:50s} {r.residual:+6.1f}   {r.same_sign:.2f}"
            f"   t {r.t_adults:+6.1f}"
        )


def step4_where(
    inputs: Inputs,
    columns: list[np.ndarray],
    splits: list[tuple[list[int], list[int]]],
) -> pd.DataFrame:
    """Which structures carry the leftover, how steadily, and which genes follow it.

    Returns residual_by_structure.csv's table.
    """
    print("\nstep 4  where the leftover lives")
    y = full_map(inputs.nano)
    leftover = residual(y, columns)
    same = same_sign_share(inputs.nano, columns, splits, leftover)
    per_adult = per_adult_leftovers(inputs.nano, columns)
    table = residual_table(inputs, y, y - leftover, leftover, same, per_adult)
    table.to_csv(RESIDUALS, index=False)
    print(
        "  above prediction (ranks; share of half-cohorts with the same sign; t across"
        " the adults' own leftovers):"
    )
    print_leftover_rows(table.head(8))
    print("  below:")
    print_leftover_rows(table.tail(6))

    # the leftover by division, with its spread over resampled cohorts
    boot = leftover_boot(inputs.nano, columns)
    divisions = division_table(inputs, leftover, boot)
    divisions.to_csv(RESIDUAL_DIVISIONS, index=False)
    division = [inputs.division.get(s, "") for s in inputs.structures]
    between = between_share(leftover, division)
    print(
        f"  by division (mean leftover in ranks; 95% over resampled adults); "
        f"{between:.0%} of its variance lies between divisions:"
    )
    for r in divisions.itertuples():
        print(
            f"    {r.division:10s} {r.n_structures:3d} structures  {r.mean:+6.1f} "
            f"({r.lo:+.1f} to {r.hi:+.1f})"
        )

    genes, sets, surr, null, tested = leftover_genes(inputs, leftover, columns, boot)
    genes.to_csv(LEFTOVER_GENES, index=False)
    sets.to_csv(LEFTOVER_SETS, index=False)
    np.savez(
        LEFTOVER_NULL,
        structures=np.array(inputs.structures),
        surrogates=surr.astype(np.float32),
        genes=np.array(tested),
        null_rho=null,
    )
    print(
        f"  every gene against the leftover ({len(genes)} genes, {len(surr)} "
        "surrogates of the leftover):"
    )
    for r in genes.head(8).itertuples():
        print(
            f"    {r.symbol:10s} rho {r.rho:+.3f}  spatial p {r.p_spatial:.4f}  "
            f"q {r.q_all:.3f}  {r.in_model}"
        )
    for g in gene_ranking.GAP_GENES:
        r = genes.set_index("symbol").loc[g]
        print(
            f"    {g}: rho {r['rho']:+.3f}, rank {int(r['rank_all'])} of {len(genes)}, "
            f"spatial p {r['p_spatial']:.4f}"
        )
    for r in sets.itertuples():
        p = f"p {r.p_spatial:.3f}" if r.tested else "not tested"
        print(
            f"    set {r.gene_set:20s} {r.n_genes:3d} genes, "
            f"median {r.median_rho:+.3f}, {p}"
        )
    return table


def main() -> dict:
    """Run the five steps on the ten adults and write their tables.

    Returns what the working figures draw: the structures used, each split's
    agreement of the two half-maps and of their leftovers, the ceiling, the
    replication the ceiling and the fit imply, the partition, and the leftover by
    structure.
    """
    OUT.mkdir(parents=True, exist_ok=True)
    inputs = load_inputs()
    step0_structures(inputs)
    print(f"  {len(inputs.structures)} structures, {len(ADULTS)} adults")

    # the ceiling, the partition, the leftover's replication, and where it lives
    splits = half_splits()
    map_agreement, explainable = step1_ceiling(inputs.nano, splits)
    covariates, _, _ = covariates_for(inputs)
    y = full_map(inputs.nano)
    table, parts, _ = step2_partition(y, covariates, inputs.terms, explainable)
    columns = model_columns(covariates, inputs.terms)
    agreement, implied = step3_replication(inputs.nano, columns, splits, map_agreement)
    residuals = step4_where(inputs, columns, splits)
    return dict(
        structures_used=inputs.structures_used,
        map_agreement=map_agreement,
        explainable=explainable,
        partition=table,
        parts=parts,
        agreement=agreement,
        implied=implied,
        residuals=residuals,
    )
