"""The adult nano map against every gene of the panel, structure by structure.

The question is whether zref can replace the old route's normalisation. The old
MATLAB route reached its adult map through a per-mouse affine fit and a
within-brain z-score; zref (volumes.cohort) is one documented transform. If the
gene ranking comes out the same, the simpler transform is enough. So the method is
the simplest one that answers it:

    nano side    one value per structure and adult, from
                 young_vs_adult/region_means_per_mouse.csv, averaged over the adults
    gene side    the region table of ish.regions
    join         on structure name, the key of both tables
    statistic    Spearman over structures, per gene and reading

Every reading is correlated, since the table holds them all: ratio and sepratio
(nano per unit autofluorescence and per unit SEP), cref, subref and zref. sepratio
was meant as the surface fraction, but the green channel is mostly
autofluorescence (adult.sep_channel_check), so it is not one; ish.arms compares
the channels.

These are correlations, not tests: two brain maps agree partly because everything
is high in cortex and hippocampus, and the usual null is inflated about 875-fold in
mouse (Fulcher 2021, docs/adult_ish_design.md). The ranking is descriptive, and the
figure says so. The old ranking (P9's frozen gene_panel_summary.csv) is compared
with the new one when it is there.

Run by run_ish_compare.py.
"""

import csv
import os
from collections import defaultdict

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

from sepmap.config import DATA

NANO = os.path.join(
    DATA, "comparisons_v2", "young_vs_adult", "region_means_per_mouse.csv"
)
GENES = os.path.join(DATA, "adult_v2", "ish", "gene_region_table.csv")
OLD = os.path.join(
    DATA,
    "comparisons",
    "merged_naive_rws_vs_ish_summary_nosmooth",
    "gene_panel_summary.csv",
)
OUT = os.path.join(DATA, "adult_v2", "ish")

# the adult groups, pooled
ADULT_GROUPS = ("naive", "rws")
READINGS = ("zref", "cref", "subref", "ratio", "sepratio")

# 200 um voxels a structure's gene value needs
MIN_ISH_VOXELS = 10

# structures a gene must share with the nano map to be correlated
MIN_STRUCTURES = 50

# genes the old and new rankings must share to be compared
MIN_GENES = 20

# the categories that carry the prediction, from gene_targets.csv
MACHINERY = ("auxiliary", "trafficking", "scaffold")


def adult_profile():
    """{reading: {structure: mean over the adults}} from the per-mouse table.

    A structure's mean is over the adults that have it, however many they are.
    """
    per = defaultdict(lambda: defaultdict(list))
    with open(NANO, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r["group"] in ADULT_GROUPS:
                per[r["reading"]][r["structure"]].append(float(r["log2_value"]))
    return {
        reading: {s: float(np.mean(v)) for s, v in d.items()}
        for reading, d in per.items()
    }


def gene_profiles():
    """{gene: {structure: expression}} and {gene: category}.

    Only structures covered by at least MIN_ISH_VOXELS voxels of the gene's grid.
    """
    out, cat = defaultdict(dict), {}
    with open(GENES, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if int(r["n_voxels"]) >= MIN_ISH_VOXELS:
                out[r["symbol"]][r["structure"]] = float(r["ish_mean"])
                cat[r["symbol"]] = r["category"]
    return out, cat


def correlate(nano, genes, category):
    """One row per reading and gene: Spearman over the structures they share.

    A gene sharing fewer than MIN_STRUCTURES structures with a reading has no row.
    """
    rows = []
    for reading in READINGS:
        profile = nano[reading]
        for gene, expr in sorted(genes.items()):
            shared = sorted(set(profile) & set(expr))
            if len(shared) < MIN_STRUCTURES:
                continue
            rho, _ = spearmanr([profile[s] for s in shared], [expr[s] for s in shared])
            rows.append(
                dict(
                    reading=reading,
                    symbol=gene,
                    category=category[gene],
                    rho=f"{rho:.4f}",
                    n_structures=len(shared),
                )
            )
    return rows


def old_ranking():
    """The MATLAB route's distance-weighted Spearman, {gene: rho}, if it is there."""
    if not os.path.exists(OLD):
        return {}
    with open(OLD, newline="", encoding="utf-8") as fh:
        return {
            r["symbol"]: float(r["r_spearman_dw"])
            for r in csv.DictReader(fh)
            if r["r_spearman_dw"] not in ("", "NaN")
        }


def rank_of(sel, gene):
    """(rank from the top, rho) of one gene in a sorted list of rows.

    (None, NaN) when the gene is not in the list.
    """
    for i, r in enumerate(sel, 1):
        if r["symbol"] == gene:
            return i, float(r["rho"])
    return None, np.nan


def report(rows, old, category):
    """Print what the decision rests on, in the order it gets asked."""
    # each reading's rows, highest rho first
    per_reading = {
        reading: sorted(
            [r for r in rows if r["reading"] == reading], key=lambda r: -float(r["rho"])
        )
        for reading in READINGS
    }

    # the top of each ranking, and where Cacng8 and Gria1 sit
    print("\ntop of each ranking, and where the two named genes sit")
    for reading, sel in per_reading.items():
        top = ", ".join(f"{r['symbol']} {float(r['rho']):+.2f}" for r in sel[:5])
        print(f"  {reading:9s} {top}")
        for gene in ("Cacng8", "Gria1"):
            i, rho = rank_of(sel, gene)
            if i:
                print(f"  {'':9s}   {gene:8s} rank {i:3d}/{len(sel)}   rho {rho:+.3f}")

    # the new ordering of the genes against the old one
    if old:
        print("\ndoes the simple normalisation reproduce the old ranking?")
        for reading, sel in per_reading.items():
            new = {r["symbol"]: float(r["rho"]) for r in sel}
            both = sorted(set(old) & set(new))
            if len(both) < MIN_GENES:
                continue
            rho, _ = spearmanr([old[g] for g in both], [new[g] for g in both])
            print(
                f"  {reading:9s} rho = {rho:+.3f} between the two orderings, "
                f"{len(both)} genes in common"
            )

    # the old route's claim, measured as it was stated: machinery genes top the list,
    # Gria1 well down it (Gria1's gap to the best machinery gene, the count above it)
    print("\nthe dissociation, stated the way it was stated of the old route")
    print(
        f"  {'reading':9s} {'best machinery':>22s} {'Gria1':>15s} {'gap':>6s}  "
        f"machinery above Gria1"
    )
    for reading, sel in per_reading.items():
        mach = [r for r in sel if r["category"] in MACHINERY]
        gria1_rank, gria1_rho = rank_of(sel, "Gria1")
        if not mach or gria1_rank is None:
            continue
        best, above = mach[0], sum(1 for r in mach if float(r["rho"]) > gria1_rho)
        print(
            f"  {reading:9s} {best['symbol']:>12s} {float(best['rho']):+8.3f} "
            f"{'rank ' + str(gria1_rank):>9s} {gria1_rho:+6.3f} "
            f"{float(best['rho']) - gria1_rho:6.3f}  "
            f"{above:>13d} / {len(mach)}"
        )

    # the same line for the old ranking
    if old:
        sel = sorted(
            (
                {"symbol": k, "rho": f"{v:.4f}", "category": category.get(k, "")}
                for k, v in old.items()
            ),
            key=lambda r: -float(r["rho"]),
        )
        gria1_rank, gria1_rho = rank_of(sel, "Gria1")
        mach = [r for r in sel if r["category"] in MACHINERY]
        above = sum(1 for r in mach if float(r["rho"]) > gria1_rho)
        print(
            f"  {'old':9s} {mach[0]['symbol']:>12s} {float(mach[0]['rho']):+8.3f} "
            f"{'rank ' + str(gria1_rank):>9s} {gria1_rho:+6.3f} "
            f"{float(mach[0]['rho']) - gria1_rho:6.3f}  "
            f"{above:>13d} / {len(mach)}"
        )


def figure(rows, old, category):
    """Draw each gene's old rho against its new one, a panel per reading.

    Machinery genes in red; saved as ish_old_vs_new.png.
    """
    # the readings in the rows, and one range for every panel
    readings = [
        r
        for r in ("zref", "cref", "ratio", "sepratio")
        if any(x["reading"] == r for x in rows)
    ]
    fig, axes = plt.subplots(
        1, len(readings), figsize=(3.1 * len(readings), 3.4), sharex=True, sharey=True
    )
    axes = np.atleast_1d(axes)

    everything = [old[g] for g in old] + [
        float(r["rho"]) for r in rows if r["reading"] in readings
    ]
    lim = [min(everything) - 0.06, max(everything) + 0.06]

    for ax, reading in zip(axes, readings):
        new = {r["symbol"]: float(r["rho"]) for r in rows if r["reading"] == reading}
        both = sorted(set(old) & set(new))
        x = np.array([old[g] for g in both])
        y = np.array([new[g] for g in both])
        mach = np.array([category.get(g, "") in MACHINERY for g in both])

        # the axes, the diagonal, the other genes, then the machinery genes
        ax.axhline(0, color="0.85", lw=0.6, zorder=0)
        ax.axvline(0, color="0.85", lw=0.6, zorder=0)
        ax.plot(lim, lim, color="0.75", lw=0.8, ls="--", zorder=1)
        ax.scatter(
            x[~mach],
            y[~mach],
            s=14,
            facecolor="0.75",
            edgecolor="0.45",
            linewidth=0.4,
            zorder=2,
        )
        ax.scatter(
            x[mach],
            y[mach],
            s=20,
            facecolor="#c0392b",
            edgecolor="0.2",
            linewidth=0.4,
            zorder=3,
        )
        for gene, dx, dy in (("Cacng8", -34, 7), ("Gria1", 6, -3)):
            if gene in both:
                ax.annotate(
                    gene,
                    (old[gene], new[gene]),
                    textcoords="offset points",
                    xytext=(dx, dy),
                    fontsize=7,
                )

        rho, _ = spearmanr(x, y)
        ax.set_title(f"{reading}\nrho = {rho:+.3f}, n = {len(both)}", fontsize=9)
        ax.set_xlabel("old route (affine + z-score)", fontsize=8)
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.tick_params(labelsize=7)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    axes[0].set_ylabel("v2 route, one reading", fontsize=8)

    fig.suptitle(
        "Gene ranking, old normalisation against v2   "
        "(red = auxiliary, trafficking or scaffold; descriptive, not tested)",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    path = os.path.join(OUT, "ish_old_vs_new.png")
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f"\n{path}")


def main():
    """Correlate every gene with the adult map, write the table, report and draw."""
    # the adult profiles and the gene profiles
    nano = adult_profile()
    genes, category = gene_profiles()
    missing = [r for r in READINGS if r not in nano]
    if missing:
        raise SystemExit(f"readings missing from {NANO}: {missing}")
    print(
        f"{len(genes)} genes, {len(nano['zref'])} adult structures, "
        f"readings {list(READINGS)}"
    )

    # one Spearman per gene and reading
    rows = correlate(nano, genes, category)
    if not rows:
        raise SystemExit("no gene shared enough structures with the nano table")

    path = os.path.join(OUT, "gene_correlations.csv")
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} rows -> {path}")

    # the report, and the comparison with the old ranking when it is there
    old = old_ranking()
    report(rows, old, category)
    if old:
        figure(rows, old, category)
