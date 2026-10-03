"""Localisation genes against expression-matched controls, on the ontology panel.

ish.roles ran the subunit-against-localisation contrast on the 100-gene panel and
found the predicted direction with no evidence for it. Its statistic compared a
4-gene composite with a 15-gene composite, so the two sides differed in noise as
well as in meaning, and more genes cannot fix that: the mouse has exactly four AMPA
subunit genes. So the question is turned round: for every gene that is not a
subunit, how much of the nano map does it explain once the subunit composite is
removed,

    partial rho(map, gene | subunit composite)

and the localisation genes are compared with the control genes on that number.
Both sides can now be large (84 localisation genes against 300 postsynaptic-density
controls), and the subunit composite enters once, as a covariate, where its noise
hurts both sides alike. If the map reports surface receptor, the genes that control
receptor localisation keep more than other postsynaptic genes once abundance is
out; if it reports total receptor, the two sets look alike.

Controls are matched on expression, because a quiet gene correlates with nothing
and reliability rises steeply with expression (ish.reliability): each localisation
gene is paired, greedily, with the unused control closest in median energy, and
the comparison with all controls is reported beside it. Labels are permuted
between the two sets, which holds set size and the expression distribution fixed.
Genes within a set are still co-expressed, so the test stays anticonservative
(Fulcher 2021), though far less than a set tested against zero, since both sets
come from the same postsynaptic population. A positive control runs the same test
on a difference that must exist (reproducible against unreproducible control
genes): if it fails, the test detects nothing and a null result means nothing.

One seeded generator feeds the permutation tests in the order they are printed, so
each p value depends on that order; the figure has a generator of its own.

Run by run_ish_panel_test.py.
"""

import csv
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import rankdata, spearmanr

from sepmap.config import DATA, SETTINGS
from sepmap.plotting import RED, tidy

# the adult groups and the reading; the structures a gene needs, the reliability of
# the sensitivity run, the permutations and the genes each check needs
ISH = SETTINGS["ish"]
ISH_PANEL_TEST = SETTINGS["ish_panel_test"]

NANO = DATA / "comparisons_v2" / "young_vs_adult" / "region_means_per_mouse.csv"
OUT = DATA / "adult_v2" / "ish"
MERGED = OUT / "gene_region_table_merged.csv"
RELIABILITY = OUT / "gene_reliability.csv"

# the adult groups, pooled, as a tuple
ADULT_GROUPS = tuple(ISH["adult_groups"])


def adult_profile(reading: str) -> dict[str, float]:
    """{structure: mean over the adults} of one reading, from the per-mouse table.

    A structure's mean is over the adults that have it, however many they are.
    """
    per = defaultdict(list)
    with open(NANO, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r["group"] in ADULT_GROUPS and r["reading"] == reading:
                per[r["structure"]].append(float(r["log2_value"]))
    return {s: float(np.mean(v)) for s, v in per.items()}


def merged_profiles() -> tuple[dict[str, dict[str, float]], dict[str, str]]:
    """{gene: {structure: mean rank}} and {gene: role}, from the merged table."""
    per, role = defaultdict(dict), {}
    with open(MERGED, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            per[r["symbol"]][r["structure"]] = float(r["rank_mean"])
            role[r["symbol"]] = r["role"]
    return per, role


def reliability_and_level() -> tuple[dict[str, float], dict[str, float]]:
    """{gene: reliability}, for the genes that have one, and {gene: median energy}."""
    rel, level = {}, {}
    with open(RELIABILITY, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r["reliability"]:
                rel[r["symbol"]] = float(r["reliability"])
            level[r["symbol"]] = float(r["median_energy"])
    return rel, level


def partial(x: np.ndarray, y: np.ndarray, z: list[np.ndarray]) -> float:
    """Spearman of `x` and `y` with the columns of `z` removed, on ranks.

    Pearson on the rank residuals; NaN when either residual is constant.
    """
    rx, ry = rankdata(x).astype(float), rankdata(y).astype(float)
    rz = np.column_stack([rankdata(c) for c in z] + [np.ones(len(rx))])

    def res(r):
        """`r` minus its least-squares fit on the ranks of `z`."""
        return r - rz @ np.linalg.lstsq(rz, r, rcond=None)[0]

    a, b = res(rx), res(ry)
    if a.std() == 0 or b.std() == 0:
        return np.nan
    return float(np.corrcoef(a, b)[0, 1])


def greedy_match(
    targets: list[str], pool: list[str], level: dict[str, float]
) -> dict[str, str]:
    """Pair each target with the closest unused control on log expression.

    Targets go in order of expression, highest first; returns {target: control}.
    """
    used, pairs = set(), {}
    order = sorted(targets, key=lambda g: -level.get(g, 0))
    for g in order:
        lg = np.log10(level.get(g, 0) + 1e-3)
        best, gap = None, np.inf
        for c in pool:
            if c in used:
                continue
            d = abs(np.log10(level.get(c, 0) + 1e-3) - lg)
            if d < gap:
                best, gap = c, d
        if best is not None:
            used.add(best)
            pairs[g] = best
    return pairs


def two_sample(
    a: np.ndarray, b: np.ndarray, rng: np.random.Generator
) -> tuple[float, float, np.ndarray]:
    """Median difference of `a` and `b`, its label-permutation p, and the null.

    The two sets are pooled and shuffled ish_panel_test.n_perm times with `rng`; p
    is two-sided,
    one added to both counts.
    """
    obs = float(np.median(a) - np.median(b))
    pool = np.concatenate([a, b])
    n = len(a)
    n_perm = ISH_PANEL_TEST["n_perm"]
    null = np.empty(n_perm)
    for i in range(n_perm):
        rng.shuffle(pool)
        null[i] = np.median(pool[:n]) - np.median(pool[n:])
    p = float((np.sum(np.abs(null) >= abs(obs)) + 1) / (n_perm + 1))
    return obs, p, null


def gene_rows(
    loc: list[str],
    ctrl: list[str],
    common: list[str],
    y: np.ndarray,
    s_comp: np.ndarray,
    expr: dict[str, dict[str, float]],
    role: dict[str, str],
    rel: dict[str, float],
    level: dict[str, float],
) -> list[dict]:
    """Plain and partial rho of every localisation and control gene, a row each.

    A gene sharing fewer than ish_panel_test.min_structures structures with the map
    and the subunit composite has no row.
    """
    rows = []
    for gene in loc + ctrl:
        shared = [s for s in common if s in expr[gene]]
        if len(shared) < ISH_PANEL_TEST["min_structures"]:
            continue
        idx = [common.index(s) for s in shared]
        g = np.array([expr[gene][s] for s in shared])
        rho, _ = spearmanr(y[idx], g)
        rows.append(
            dict(
                symbol=gene,
                role=role[gene],
                n_structures=len(shared),
                rho=float(rho),
                rho_partial=partial(y[idx], g, [s_comp[idx]]),
                reliability=rel.get(gene, np.nan),
                median_energy=level.get(gene, np.nan),
            )
        )
    return rows


def write_table(rows: list[dict]) -> None:
    """Write panel_test.csv, the per-gene table, floats to four decimals."""
    path = OUT / "panel_test.csv"
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=[
                "symbol",
                "role",
                "n_structures",
                "rho",
                "rho_partial",
                "reliability",
                "median_energy",
                "matched_to",
            ],
        )
        w.writeheader()
        w.writerows(
            {k: (f"{v:.4f}" if isinstance(v, float) else v) for k, v in r.items()}
            for r in rows
        )
    print(f"{len(rows)} genes -> {path}")


def localisation_tests(
    loc: list[str],
    ctrl: list[str],
    matched_ctrl: list[str],
    by: dict[str, dict],
    rng: np.random.Generator,
) -> dict[str, tuple]:
    """The test against all controls and against the matched ones, then plain rho.

    Returns {label: (a, b, difference, p, null)} of the two partial-rho tests;
    `rng` is drawn from in the order the tests are printed.
    """
    print("\nTEST -- partial rho given the subunit composite, localisation vs control")
    results = {}
    for label, cs in (("all controls", ctrl), ("expression-matched", matched_ctrl)):
        a = np.array([by[g]["rho_partial"] for g in loc])
        b = np.array([by[g]["rho_partial"] for g in cs])
        a, b = a[np.isfinite(a)], b[np.isfinite(b)]
        obs, p, null = two_sample(a, b, rng)
        results[label] = (a, b, obs, p, null)
        detectable = float(np.percentile(np.abs(null), 95))
        print(
            f"  {label:20s} localisation {np.median(a):+.3f} (n={len(a)}), "
            f"control {np.median(b):+.3f} (n={len(b)}), "
            f"difference {obs:+.3f}, p = {p:.4f}"
        )
        print(
            f"  {'':20s} a difference of {detectable:+.3f} would have been detected "
            f"at p < 0.05; observed is {abs(obs) / detectable:.2f} of that"
        )

    # and the plain correlation, for comparison with everything reported before
    a = np.array([by[g]["rho"] for g in loc])
    b = np.array([by[g]["rho"] for g in matched_ctrl])
    obs, p, _ = two_sample(a, b, rng)
    print(
        f"  {'before partialling':20s} localisation {np.median(a):+.3f}, "
        f"control {np.median(b):+.3f}, difference {obs:+.3f}, p = {p:.4f}"
    )
    return results


def sensitivity_test(
    loc: list[str], matched_ctrl: list[str], by: dict[str, dict], rng: np.random.Generator
) -> None:
    """Print the matched test again on the genes whose map is reproducible."""
    min_rel = ISH_PANEL_TEST["min_reliability"]
    min_genes = ISH_PANEL_TEST["min_sensitivity_genes"]
    good_loc = [g for g in loc if by[g]["reliability"] >= min_rel]
    good_ctrl = [g for g in matched_ctrl if by[g]["reliability"] >= min_rel]
    if len(good_loc) >= min_genes and len(good_ctrl) >= min_genes:
        a = np.array([by[g]["rho_partial"] for g in good_loc])
        b = np.array([by[g]["rho_partial"] for g in good_ctrl])
        obs2, p2, _ = two_sample(a, b, rng)
        print(
            f"\n  sensitivity, reliability >= {min_rel}: "
            f"{len(good_loc)} vs {len(good_ctrl)} genes, "
            f"difference {obs2:+.3f}, p = {p2:.4f}"
        )


def positive_control(
    ctrl: list[str], by: dict[str, dict], rng: np.random.Generator
) -> None:
    """Print the same test on a contrast that must exist.

    Among the control genes, those with a reproducible map should correlate better
    with anything than the others; if not, the test detects nothing at all.
    """
    have = [g for g in ctrl if np.isfinite(by[g]["reliability"])]
    if len(have) >= ISH_PANEL_TEST["min_control_genes"]:
        cut = float(np.median([by[g]["reliability"] for g in have]))
        hi = np.array([abs(by[g]["rho"]) for g in have if by[g]["reliability"] >= cut])
        lo = np.array([abs(by[g]["rho"]) for g in have if by[g]["reliability"] < cut])
        obs_c, p_c, _ = two_sample(hi, lo, rng)
        print(
            f"\n  positive control -- reproducible vs unreproducible control genes "
            f"(reliability split at {cut:.2f}):"
        )
        print(
            f"    |rho| {np.median(hi):.3f} (n={len(hi)}) against {np.median(lo):.3f} "
            f"(n={len(lo)}), difference {obs_c:+.3f}, p = {p_c:.4f}"
        )


def print_top(loc: list[str], by: dict[str, dict]) -> None:
    """Print the ten localisation genes with the highest partial rho."""
    print("\n  top localisation genes by partial rho:")
    for r in sorted((by[g] for g in loc), key=lambda r: -r["rho_partial"])[:10]:
        shown = "n/a" if np.isnan(r["reliability"]) else f"{r['reliability']:.2f}"
        print(
            f"    {r['symbol']:10s} partial {r['rho_partial']:+.3f}  "
            f"plain {r['rho']:+.3f}  reliability {shown}"
        )


def main() -> None:
    """Run the test, its sensitivity run and the positive control; write and draw."""
    # the map, the merged gene profiles, each gene's reliability and level
    nano = adult_profile(ISH["reading"])
    expr, role = merged_profiles()
    rel, level = reliability_and_level()

    # the three sets of genes
    subunit = sorted(g for g in expr if role[g] == "subunit")
    loc = sorted(g for g in expr if role[g] == "localisation")
    ctrl = sorted(g for g in expr if role[g] == "control_psd")
    print(
        f"{len(expr)} genes: {len(subunit)} subunit, {len(loc)} localisation, "
        f"{len(ctrl)} control; reading {ISH['reading']}"
    )
    print(f"  subunits: {', '.join(subunit)}")

    # one structure set for everything, so no comparison moves the ground
    common = sorted(set(nano).intersection(*[set(expr[g]) for g in subunit]))
    y = np.array([nano[s] for s in common])
    s_comp = np.mean([rankdata([expr[g][s] for s in common]) for g in subunit], axis=0)
    print(f"  {len(common)} structures shared by the map and all subunit genes")

    # plain and partial rho of every localisation and control gene
    rows = gene_rows(loc, ctrl, common, y, s_comp, expr, role, rel, level)
    by = {r["symbol"]: r for r in rows}
    loc = [g for g in loc if g in by]
    ctrl = [g for g in ctrl if g in by]
    print(f"  usable: {len(loc)} localisation, {len(ctrl)} control")

    # each localisation gene's expression-matched control, and the per-gene table
    pairs = greedy_match(loc, ctrl, level)
    for r in rows:
        r["matched_to"] = pairs.get(r["symbol"], "")
    matched_ctrl = sorted(set(pairs.values()))
    write_table(rows)

    # the tests; this generator feeds every permutation test, in this order: all
    # controls, the matched ones, before partialling, sensitivity, positive control
    rng = np.random.default_rng(0)
    results = localisation_tests(loc, ctrl, matched_ctrl, by, rng)
    sensitivity_test(loc, matched_ctrl, by, rng)
    positive_control(ctrl, by, rng)

    # the top localisation genes, and the figure
    print_top(loc, by)
    figure(results, by, loc, matched_ctrl, ctrl, level)


def panel_matched(
    ax: plt.Axes,
    a: np.ndarray,
    b: np.ndarray,
    obs: float,
    p: float,
    rng: np.random.Generator,
) -> None:
    """Draw the matched controls against the localisation genes, the test."""
    for i, (v, colour) in enumerate(((b, "0.65"), (a, RED))):
        ax.scatter(
            np.full(len(v), i) + rng.uniform(-0.14, 0.14, len(v)),
            v,
            s=14,
            facecolor=colour,
            edgecolor="0.25",
            linewidth=0.4,
            zorder=2,
        )
        ax.plot([i - 0.3, i + 0.3], [np.median(v)] * 2, color="0.15", lw=1.8, zorder=3)
    ax.axhline(0, color="0.85", lw=0.7, zorder=0)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(
        [f"matched\ncontrols ({len(b)})", f"localisation\n({len(a)})"], fontsize=8
    )
    ax.set_ylabel("partial rho with the map, subunits removed", fontsize=8)
    ax.set_title(f"the test\ndifference {obs:+.3f}, p = {p:.4f}", fontsize=9)


def panel_expression(
    ax: plt.Axes,
    ctrl: list[str],
    matched_ctrl: list[str],
    loc: list[str],
    level: dict[str, float],
    rng: np.random.Generator,
) -> None:
    """Draw the expression level of each set: why matching was needed."""
    for i, (genes, colour, label) in enumerate(
        (
            (ctrl, "0.8", "all controls"),
            (matched_ctrl, "0.5", "matched"),
            (loc, RED, "localisation"),
        )
    ):
        v = np.log10([level.get(g, 0) + 1e-3 for g in genes])
        ax.scatter(
            np.full(len(v), i) + rng.uniform(-0.14, 0.14, len(v)),
            v,
            s=10,
            facecolor=colour,
            edgecolor="0.3",
            linewidth=0.3,
            zorder=2,
        )
        ax.plot([i - 0.3, i + 0.3], [np.median(v)] * 2, color="0.15", lw=1.7, zorder=3)
    ax.set_xticks(range(3))
    ax.set_xticklabels(["all\ncontrols", "matched\ncontrols", "localisation"], fontsize=8)
    ax.set_ylabel("log10 median expression energy", fontsize=8)
    ax.set_title(
        "why matching was needed\na quiet gene correlates with nothing", fontsize=9
    )


def panel_null(ax: plt.Axes, null: np.ndarray, obs: float, p: float) -> None:
    """Draw the permutation null and the observed difference."""
    ax.hist(null, bins=60, color="0.72", edgecolor="0.35", linewidth=0.3)
    ax.axvline(obs, color=RED, lw=2)
    ax.set_xlabel("median(localisation) - median(control)", fontsize=8)
    ax.set_ylabel(f"label permutations ({len(null):,})", fontsize=8)
    ax.set_title(f"the null\np = {p:.4f}", fontsize=9)


def figure(
    results: dict[str, tuple],
    by: dict[str, dict],
    loc: list[str],
    matched_ctrl: list[str],
    ctrl: list[str],
    level: dict[str, float],
) -> None:
    """Draw the matched test, the expression of each set and the null.

    Saved as ish_panel_test.png; the jitter has a generator of its own.
    """
    fig, axes = plt.subplots(1, 3, figsize=(13.8, 4.3))
    rng = np.random.default_rng(0)

    # left: the matched controls against the localisation genes; middle: the
    # expression level of each set; right: the permutation null and the observed
    # difference
    a, b, obs, p, null = results["expression-matched"]
    panel_matched(axes[0], a, b, obs, p, rng)
    panel_expression(axes[1], ctrl, matched_ctrl, loc, level, rng)
    panel_null(axes[2], null, obs, p)

    for ax in axes:
        tidy(ax)

    # the title states the question
    fig.suptitle(
        "Once receptor abundance is taken out of the map, do the genes that "
        "control receptor localisation explain more than other postsynaptic "
        "genes do?",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    path = OUT / "ish_panel_test.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f"\n{path}")
