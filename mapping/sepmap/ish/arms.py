"""The channel arms against the genes: membrane pool or total receptor?

The gene ranking alone cannot tell the two apart: Gria1, Cacng8, Dlg2 and Grip1 are
all postsynaptic, so no grouping of genes separates the surface pool from total
receptor, and both predict a postsynaptic map. The distinction was to come from
the channels, which adult.arms puts on one footing:

    sepauto     SEP / auto     meant as total receptor, wherever it sits
    ratio       nano / auto    receptor at the surface
    sepratio    nano / SEP     meant as the surface fraction, auto divided out

Two tests, both stated before the numbers were computed. Test 1, the order across
arms: if SEP reports total receptor and nano the surface pool, then against Gria1
the arms order sepauto > ratio > sepratio, and against trafficking and anchoring
machinery the order reverses; reported per gene and as sepratio - sepauto. Test 2,
the partial correlation given Gria1: even a perfect surface map correlates with
Gria1, since there is no surface receptor without receptor, so the question is
what is left once Gria1 is taken out, rho(arm, gene | Gria1) on ranks. If the map
reports abundance only, nothing survives; if it reports the surface pool, the
machinery genes still predict the residual.

The premise of test 1 failed: SEP / auto correlates with Gria1 at -0.11, and
adult.sep_channel_check shows that in this fixed, cleared tissue the green channel
is mostly autofluorescence (rho 0.79 +- 0.04 with the autofluorescence channel
across the ten adults, against 0.26 for nano). So sepauto is not total receptor
and sepratio is not a surface fraction (it follows ratio at rho 0.89 to 0.97 in
every mouse): test 1 measures nothing about membrane against total, and stays as
the record of how that was found. Test 2 holds, since it does not use the green
channel. Gria1 expression stands in for total receptor, a weaker reference than a
total-protein channel: mRNA is not protein, and a residual that trafficking genes
predict also fits a regional translation gradient.

Gria1's map is one ISH experiment, and a noisy one would weaken every correlation
with it and inflate what survives the partial, so its coverage is printed. The p
values that would go with these rho are inflated (Fulcher 2021), so the arms are
compared with each other on the same genes and structures, a paired contrast in
which that inflation largely cancels, never against zero.

Run by run_ish_arms.py.
"""

import csv
import os
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import rankdata, spearmanr, wilcoxon

from sepmap.config import DATA, SETTINGS
from sepmap.ish.compare import gene_profiles
from sepmap.plotting import DARK_BLUE, RED, tidy

# the gene that stands in for total receptor, partialled out in test 2, and the
# structures a gene must share with the arms (and with that gene)
ISH = SETTINGS["ish"]

ARMS_CSV = os.path.join(DATA, "adult_v2", "arms", "region_means_arms.csv")
OUT = os.path.join(DATA, "adult_v2", "arms")

# the arms in the order of test 1, and their labels in the figure
ARMS = ("sepauto", "ratio", "sepratio")
LABEL = {
    "sepauto": "SEP / auto\ntotal receptor",
    "ratio": "nano / auto\nsurface receptor",
    "sepratio": "nano / SEP\nsurface fraction",
}

# the genes the prediction is about: AMPAR anchoring and trafficking, named in
# gene_targets.csv, not chosen after seeing this result
MACHINERY = ("auxiliary", "trafficking", "scaffold")


def arm_profiles() -> dict[str, dict[str, float]]:
    """{arm: {structure: mean over the adults}}, from region_means_arms.csv."""
    per = defaultdict(lambda: defaultdict(list))
    with open(ARMS_CSV, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            per[r["arm"]][r["structure"]].append(float(r["log2_value"]))
    return {arm: {s: float(np.mean(v)) for s, v in d.items()} for arm, d in per.items()}


def partial_spearman(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> float:
    """Spearman of `x` and `y` with `z` removed: Pearson on the rank residuals.

    NaN when either residual is constant.
    """
    rx, ry, rz = (rankdata(v).astype(float) for v in (x, y, z))
    zc = np.column_stack([rz, np.ones_like(rz)])

    def resid(r):
        """`r` minus its least-squares fit on the control gene's ranks."""
        return r - zc @ np.linalg.lstsq(zc, r, rcond=None)[0]

    a, b = resid(rx), resid(ry)
    if a.std() == 0 or b.std() == 0:
        return np.nan
    return float(np.corrcoef(a, b)[0, 1])


def correlate(
    arms: dict[str, dict[str, float]],
    genes: dict[str, dict[str, float]],
    category: dict[str, str],
    control: str,
) -> list[dict]:
    """One row per arm and gene: rho, and rho with the control gene removed.

    Every arm is restricted to the same structures, those measured in all three
    arms, in this gene and in the control gene, so the three numbers of a gene
    are a paired comparison, not three different samples. The control gene's own
    partial rho is NaN.
    """
    shared_arms = set.intersection(*[set(arms[a]) for a in ARMS])
    rows = []
    for gene, expr in sorted(genes.items()):
        if gene == control:
            common = sorted(shared_arms & set(expr))
        else:
            common = sorted(shared_arms & set(expr) & set(genes[control]))
        if len(common) < ISH["min_structures"]:
            continue

        # the gene, the control gene, and each arm, on the same structures
        y = np.array([expr[s] for s in common])
        c = np.array([genes[control][s] for s in common])
        for arm in ARMS:
            x = np.array([arms[arm][s] for s in common])
            rho, _ = spearmanr(x, y)
            rows.append(
                dict(
                    arm=arm,
                    symbol=gene,
                    category=category[gene],
                    n_structures=len(common),
                    rho=float(rho),
                    rho_partial=(
                        np.nan if gene == control else partial_spearman(x, y, c)
                    ),
                )
            )
    return rows


def by_gene(rows: list[dict], field: str) -> dict[str, dict[str, float]]:
    """{gene: {arm: value}} for one column."""
    out = defaultdict(dict)
    for r in rows:
        out[r["symbol"]][r["arm"]] = r[field]
    return out


def report(
    rows: list[dict], category: dict[str, str]
) -> tuple[dict[str, dict[str, float]], dict[str, dict[str, float]], list[str]]:
    """Print both tests and the paired contrast; return the rho tables and machinery.

    Returns ({gene: {arm: rho}}, {gene: {arm: partial rho}}, the machinery genes
    sorted by their sepratio rho, highest first).
    """
    control = ISH["control_gene"]
    plain, partial = by_gene(rows, "rho"), by_gene(rows, "rho_partial")
    mach = sorted(
        (g for g in plain if category[g] in MACHINERY),
        key=lambda g: -plain[g]["sepratio"],
    )

    # test 1: Gria1 and the top machinery genes across the arms
    print("\nTEST 1 -- the ordering across arms")
    print(
        f"  {'':12s} {'sepauto':>9s} {'ratio':>9s} {'sepratio':>9s}   "
        f"{'sepratio-sepauto':>17s}"
    )
    named = [control] + mach[:6]
    for gene in named:
        rhos = plain[gene]
        print(
            f"  {gene:12s} {rhos['sepauto']:+9.3f} {rhos['ratio']:+9.3f} "
            f"{rhos['sepratio']:+9.3f}   {rhos['sepratio'] - rhos['sepauto']:+17.3f}"
        )

    # test 1: the swing towards the surface fraction, Gria1 against the machinery
    swing = {g: plain[g]["sepratio"] - plain[g]["sepauto"] for g in plain}
    mach_swing = [swing[g] for g in plain if category[g] in MACHINERY]
    print(f"\n  {control} swing towards the surface fraction: {swing[control]:+.3f}")
    print(
        f"  machinery genes (n = {len(mach_swing)}):  "
        f"median {np.median(mach_swing):+.3f}, "
        f"{sum(1 for v in mach_swing if v > swing[control])} of {len(mach_swing)} "
        f"above {control}"
    )
    print(
        "  the prediction is that Gria1 swings DOWN and the machinery swings up "
        "or falls less"
    )

    # test 2: what is left of the machinery genes with Gria1 partialled out
    print(f"\nTEST 2 -- with {control} partialled out, what is left")
    print(f"  {'':12s} {'sepauto':>9s} {'ratio':>9s} {'sepratio':>9s}")
    for gene in mach[:8]:
        rhos = partial[gene]
        print(
            f"  {gene:12s} {rhos['sepauto']:+9.3f} {rhos['ratio']:+9.3f} "
            f"{rhos['sepratio']:+9.3f}"
        )
    for arm in ARMS:
        vals = [
            partial[g][arm]
            for g in plain
            if category[g] in MACHINERY and np.isfinite(partial[g][arm])
        ]
        print(
            f"  {arm:12s} machinery median {np.median(vals):+.3f}, "
            f"{sum(1 for v in vals if v > 0)} of {len(vals)} positive"
        )

    # the paired contrast: same genes, same structures, two arms
    a = np.array([plain[g]["sepratio"] for g in mach])
    b = np.array([plain[g]["sepauto"] for g in mach])
    stat, p = wilcoxon(a, b)
    print(
        f"\n  machinery, sepratio vs sepauto, paired over {len(mach)} genes: "
        f"median difference {np.median(a - b):+.3f}, Wilcoxon p = {p:.4f}"
    )
    return plain, partial, mach


def figure(
    plain: dict[str, dict[str, float]],
    partial: dict[str, dict[str, float]],
    mach: list[str],
    category: dict[str, str],
) -> None:
    """Draw the two tests in three panels; saved as arms_vs_genes.png."""
    control = ISH["control_gene"]
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.3))

    # left: every gene's rho across the three arms, Gria1 and machinery picked out
    ax = axes[0]
    x = np.arange(len(ARMS))
    for gene in sorted(plain):
        if gene == control or category[gene] in MACHINERY:
            continue
        ax.plot(x, [plain[gene][a] for a in ARMS], color="0.85", lw=0.7, zorder=1)
    for gene in mach:
        ax.plot(
            x,
            [plain[gene][a] for a in ARMS],
            color=RED,
            lw=0.9,
            alpha=0.55,
            zorder=2,
        )
    ax.plot(
        x,
        [plain[control][a] for a in ARMS],
        color=DARK_BLUE,
        lw=2.4,
        zorder=4,
        marker="o",
        ms=5,
    )
    ax.annotate(
        control,
        (2, plain[control]["sepratio"]),
        color=DARK_BLUE,
        fontsize=9,
        xytext=(6, -2),
        textcoords="offset points",
    )
    ax.set_xticks(x)
    ax.set_xticklabels([LABEL[a] for a in ARMS], fontsize=7.5)
    ax.set_ylabel("Spearman with the arm, over structures", fontsize=8)
    ax.set_title(
        "Test 1: every gene across the three arms\n"
        "red = auxiliary / trafficking / scaffold",
        fontsize=9,
    )

    # middle: the swing, machinery against everything else
    ax = axes[1]
    swing = {g: plain[g]["sepratio"] - plain[g]["sepauto"] for g in plain}
    groups = [
        [swing[g] for g in plain if category[g] not in MACHINERY and g != control],
        [swing[g] for g in mach],
    ]
    rng = np.random.default_rng(0)
    for i, (vals, colour) in enumerate(zip(groups, ("0.65", RED))):
        ax.scatter(
            np.full(len(vals), i) + rng.uniform(-0.12, 0.12, len(vals)),
            vals,
            s=15,
            facecolor=colour,
            edgecolor="0.25",
            linewidth=0.4,
            zorder=2,
        )
        ax.plot(
            [i - 0.28, i + 0.28], [np.median(vals)] * 2, color="0.15", lw=1.7, zorder=3
        )
    ax.axhline(swing[control], color=DARK_BLUE, lw=1.4, ls="--", zorder=1)
    ax.annotate(
        control,
        (1.35, swing[control]),
        color=DARK_BLUE,
        fontsize=8,
        va="bottom",
        ha="right",
    )
    ax.axhline(0, color="0.8", lw=0.7, zorder=0)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(
        [f"other genes\n({len(groups[0])})", f"machinery\n({len(groups[1])})"], fontsize=8
    )
    ax.set_ylabel("rho(surface fraction) - rho(total receptor)", fontsize=8)
    ax.set_title("Test 1: the swing towards the surface fraction", fontsize=9)

    # right: what survives partialling out Gria1
    ax = axes[2]
    for i, arm in enumerate(ARMS):
        vals = [partial[g][arm] for g in mach if np.isfinite(partial[g][arm])]
        ax.scatter(
            np.full(len(vals), i) + rng.uniform(-0.12, 0.12, len(vals)),
            vals,
            s=15,
            facecolor=RED,
            edgecolor="0.25",
            linewidth=0.4,
            zorder=2,
        )
        ax.plot(
            [i - 0.28, i + 0.28], [np.median(vals)] * 2, color="0.15", lw=1.7, zorder=3
        )
    ax.axhline(0, color="0.8", lw=0.7, zorder=0)
    ax.set_xticks(range(len(ARMS)))
    ax.set_xticklabels([LABEL[a] for a in ARMS], fontsize=7.5)
    ax.set_ylabel(f"partial rho with {control} removed", fontsize=8)
    ax.set_title(f"Test 2: machinery genes, {control} partialled out", fontsize=9)

    for ax in axes:
        tidy(ax)

    # the title says how the panels compare
    fig.suptitle(
        "Membrane pool or total receptor?  The contrast is between channels, "
        "not between genes.  Descriptive: arms are compared with each other on the "
        "same genes and structures, never against zero.",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    path = os.path.join(OUT, "arms_vs_genes.png")
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f"\n{path}")


def main() -> None:
    """Correlate every gene with every arm, write the table, report and draw."""
    control = ISH["control_gene"]
    # the arm profiles and the gene profiles, the control gene among them
    arms = arm_profiles()
    missing = [a for a in ARMS if a not in arms]
    if missing:
        raise ValueError(
            f"arms missing from {ARMS_CSV}: {missing}. Run run_adult_arms.py"
        )
    genes, category = gene_profiles()
    if control not in genes:
        raise ValueError(f"{control} is not in the gene table; it is the control here")
    print(
        f"{len(genes)} genes; {control} measured in {len(genes[control])} structures; "
        f"arms {ARMS}"
    )

    # plain and partial rho per arm and gene
    rows = correlate(arms, genes, category, control)
    path = os.path.join(OUT, "arm_gene_correlations.csv")
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(
            {k: (f"{v:.4f}" if isinstance(v, float) else v) for k, v in r.items()}
            for r in rows
        )
    print(f"{len(rows)} rows -> {path}")

    # the two tests, printed and drawn
    plain, partial, mach = report(rows, category)
    figure(plain, partial, mach, category)
