"""How reliable one Allen ISH map is, measured on the genes measured twice.

Every correlation here is between the nano map and one ISH experiment per gene:
one mouse, one set of sections, one staining run. A correlation is capped by the
reliability of both things correlated, so a gene with a noisy map looks
uninformative whether it is or not. Many genes of the ontology panel have a
coronal and a sagittal experiment, two independent passes at the same gene in
different animals and sectioning, and the correlation of the two profiles is a
test-retest reliability of one Allen map. Two things come of it:

    reliability    per gene, the median Spearman between its experiments: the cap
                   on every correlation with that gene, and a quality filter that
                   is measured rather than assumed
    merged map     the mean of a gene's experiments, a better estimate than either,
                   averaged on ranks because expression energy is not calibrated
                   across experiments (two runs of one gene can differ by a factor
                   that has nothing to do with the brain)

Coronal against sagittal is the fair pairing and also the harshest: the two planes
are reconstructed differently and a sagittal experiment covers one hemisphere, so
their agreement is a lower bound on what two experiments of one plane would give.
Pairs of the same plane are reported apart, and the difference between the two
kinds is worth a look before leaning on either.

Run by run_ish_reliability.py.
"""

import csv
import itertools
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import rankdata, spearmanr

from sepmap.config import DATA, SETTINGS
from sepmap.plotting import RED, tidy

OUT = DATA / "adult_v2" / "ish"
MERGED = OUT / "gene_region_table_merged.csv"

# the table to read is the one run_ish_regions wrote for a panel pass of
# settings.toml ([ish_panels]); by default the ontology panel, whose genes include
# those with more than one experiment
ISH_PANELS = SETTINGS["ish_panels"]
DEFAULT_PANEL = "ontology"

# the voxels a gene value needs, and the structures two experiments must share
ISH = SETTINGS["ish"]


def load(
    path: Path,
) -> tuple[
    dict[str, dict[str, dict[str, float]]], dict[tuple[str, str], tuple[str, str]]
]:
    """Each experiment's profile, and each experiment's plane and role.

    Returns ({gene: {experiment: {structure: energy}}}, {(gene, experiment):
    (plane, role)}), structures with fewer than ish.min_voxels voxels left out.
    """
    per = defaultdict(lambda: defaultdict(dict))
    meta = {}
    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if int(r["n_voxels"]) < ISH["min_voxels"]:
                continue
            per[r["symbol"]][r["experiment_id"]][r["structure"]] = float(r["ish_mean"])
            meta[(r["symbol"], r["experiment_id"])] = (r.get("plane", ""), r["category"])
    return per, meta


def pair_reliability(
    profiles: dict[str, dict[str, float]],
) -> list[tuple[str, str, float, int]]:
    """(experiment a, experiment b, Spearman, n structures) for each pair of a gene.

    A pair sharing fewer than ish.min_structures_pair structures is left out.
    """
    out = []
    for a, b in itertools.combinations(sorted(profiles), 2):
        shared = sorted(set(profiles[a]) & set(profiles[b]))
        if len(shared) < ISH["min_structures_pair"]:
            continue
        rho, _ = spearmanr(
            [profiles[a][s] for s in shared], [profiles[b][s] for s in shared]
        )
        out.append((a, b, float(rho), len(shared)))
    return out


def merge(
    profiles: dict[str, dict[str, float]],
) -> tuple[dict[str, float], dict[str, int]]:
    """One profile per gene: the mean of its experiments' ranks.

    Ranks rather than values because expression energy carries an arbitrary
    per-experiment scale. Structures are kept if at least one experiment has
    them, and each experiment's ranks are put on 0-1 first so that experiments
    covering different numbers of structures can still be averaged. Returns
    ({structure: mean rank}, {structure: n experiments}).
    """
    acc = defaultdict(list)
    for exp, prof in profiles.items():
        keys = sorted(prof)
        r = rankdata([prof[k] for k in keys])
        r = (r - 1) / max(len(r) - 1, 1)
        for k, v in zip(keys, r):
            acc[k].append(v)
    return (
        {k: float(np.mean(v)) for k, v in acc.items()},
        {k: len(v) for k, v in acc.items()},
    )


def gene_tables(
    per: dict[str, dict[str, dict[str, float]]],
    meta: dict[tuple[str, str], tuple[str, str]],
) -> tuple[list[dict], list[dict], dict[str, list[float]]]:
    """Each gene's reliability, by pairing too, and its merged profile.

    Returns the rows of gene_reliability.csv and of gene_region_table_merged.csv,
    and {pairing: [rho, ...]}.
    """
    rows, merged_rows = [], []
    rel_by_pairing = defaultdict(list)
    for gene, profiles in sorted(per.items()):
        pairs = pair_reliability(profiles)
        rel = float(np.median([p[2] for p in pairs])) if pairs else np.nan
        planes = sorted({meta[(gene, e)][0] for e in profiles})
        role = meta[(gene, sorted(profiles)[0])][1]
        for a, b, rho, n in pairs:
            pa, pb = meta[(gene, a)][0], meta[(gene, b)][0]
            rel_by_pairing["same plane" if pa == pb else "coronal vs sagittal"].append(
                rho
            )

        prof, n_exp = merge(profiles)
        level = float(np.median([v for p in profiles.values() for v in p.values()]))
        rows.append(
            dict(
                symbol=gene,
                role=role,
                n_experiments=len(profiles),
                planes=" ".join(planes),
                n_pairs=len(pairs),
                reliability=("" if np.isnan(rel) else f"{rel:.4f}"),
                n_structures=len(prof),
                median_energy=f"{level:.4g}",
            )
        )
        for structure, value in prof.items():
            merged_rows.append(
                dict(
                    symbol=gene,
                    role=role,
                    structure=structure,
                    rank_mean=f"{value:.6f}",
                    n_experiments=n_exp[structure],
                )
            )
    return rows, merged_rows, rel_by_pairing


def merged_profiles() -> tuple[dict[str, dict[str, float]], dict[str, str]]:
    """{gene: {structure: mean rank}} and {gene: role}, from the merged table."""
    per, role = defaultdict(dict), {}
    with open(MERGED, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            per[r["symbol"]][r["structure"]] = float(r["rank_mean"])
            role[r["symbol"]] = r["role"]
    return per, role


def write_tables(rows: list[dict], merged_rows: list[dict]) -> None:
    """Write gene_reliability.csv and gene_region_table_merged.csv."""
    p1 = OUT / "gene_reliability.csv"
    with open(p1, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    p2 = MERGED
    with open(p2, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(merged_rows[0].keys()))
        w.writeheader()
        w.writerows(merged_rows)
    print(f"{len(rows)} genes -> {p1}\n{len(merged_rows):,} rows -> {p2}")


def print_reliability(
    rows: list[dict], rel_by_pairing: dict[str, list[float]]
) -> np.ndarray:
    """Print the distribution, by pairing, and the genes the panel test turns on.

    Returns the reliability of every gene measured more than once.
    """
    have = [r for r in rows if r["reliability"] != ""]
    rel = np.array([float(r["reliability"]) for r in have])
    print(f"\n{len(have)} genes measured more than once")
    print(
        f"  reliability  median {np.median(rel):+.3f}   "
        f"quartiles [{np.percentile(rel, 25):+.3f} {np.percentile(rel, 75):+.3f}]   "
        f"{(rel < 0.3).sum()} below 0.3"
    )
    for pairing, v in sorted(rel_by_pairing.items()):
        print(f"  {pairing:20s} n={len(v):4d}  median {np.median(v):+.3f}")

    print("\n  the genes the test turns on:")
    for gene in ("Gria1", "Gria2", "Gria3", "Gria4", "Cacng8", "Dlg4", "Dlg2", "Nsf"):
        r = next((x for x in rows if x["symbol"] == gene), None)
        if r:
            print(
                f"    {gene:8s} {r['n_experiments']} experiments ({r['planes']}), "
                f"reliability {r['reliability'] or 'n/a'}"
            )
    return rel


def main(panel_name: str = DEFAULT_PANEL) -> None:
    """Write each gene's reliability and merged profile for one panel pass; draw."""
    # the region table of the panel pass
    table_name = ISH_PANELS[panel_name]["table"]
    path = OUT / table_name
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found -- run run_ish_regions.py with --panel {panel_name} first"
        )
    per, meta = load(path)
    print(
        f"{len(per)} genes, {sum(len(v) for v in per.values())} experiments, "
        f"from {table_name}"
    )

    # each gene's reliability, by pairing too, and its merged profile, as two tables
    rows, merged_rows, rel_by_pairing = gene_tables(per, meta)
    write_tables(rows, merged_rows)

    # the distribution, by pairing, and the genes the panel test turns on
    rel = print_reliability(rows, rel_by_pairing)
    figure(rows, rel, rel_by_pairing, per, meta)


def panel_distribution(ax: plt.Axes, rel: np.ndarray) -> None:
    """Draw the distribution of the reliability."""
    ax.hist(rel, bins=40, color="0.7", edgecolor="0.35", linewidth=0.4)
    ax.axvline(float(np.median(rel)), color=RED, lw=1.6)
    ax.set_xlabel("Spearman between two experiments of the same gene", fontsize=8)
    ax.set_ylabel("genes", fontsize=8)
    ax.set_title(
        f"how reliable one Allen map is\nmedian {np.median(rel):+.2f}, n = {len(rel)}",
        fontsize=9,
    )


def panel_expression(ax: plt.Axes, rows: list[dict]) -> None:
    """Draw the reliability against the expression level."""
    have = [r for r in rows if r["reliability"] != ""]
    x = np.log10([float(r["median_energy"]) + 1e-3 for r in have])
    y = [float(r["reliability"]) for r in have]
    ax.scatter(x, y, s=10, facecolor="0.55", edgecolor="0.25", linewidth=0.3)
    ax.axhline(0, color="0.85", lw=0.7)
    ax.set_xlabel("log10 median expression energy", fontsize=8)
    ax.set_ylabel("reliability", fontsize=8)
    ax.set_title(
        f"a quiet gene is an unreliable one\nrho = {spearmanr(x, y).statistic:+.2f}",
        fontsize=9,
    )


def panel_pairing(ax: plt.Axes, rel_by_pairing: dict[str, list[float]]) -> None:
    """Draw the reliability by the pairing available; the jitter has its own generator."""
    rng = np.random.default_rng(0)
    keys = sorted(rel_by_pairing)
    for i, k in enumerate(keys):
        v = rel_by_pairing[k]
        ax.scatter(
            np.full(len(v), i) + rng.uniform(-0.12, 0.12, len(v)),
            v,
            s=11,
            facecolor="0.6",
            edgecolor="0.25",
            linewidth=0.3,
            zorder=2,
        )
        ax.plot([i - 0.28, i + 0.28], [np.median(v)] * 2, color=RED, lw=1.8, zorder=3)
    ax.set_xticks(range(len(keys)))
    ax.set_xticklabels([f"{k}\n({len(rel_by_pairing[k])})" for k in keys], fontsize=8)
    ax.axhline(0, color="0.85", lw=0.7)
    ax.set_ylabel("reliability", fontsize=8)
    ax.set_title("which pairing was available", fontsize=9)


def figure(
    rows: list[dict],
    rel: np.ndarray,
    rel_by_pairing: dict[str, list[float]],
    per: dict,
    meta: dict,
) -> None:
    """Draw the reliability, its dependence on expression, and the pairings.

    Saved as ish_reliability.png.
    """
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.0))

    # left: the distribution of the reliability; middle: reliability against
    # expression level; right: reliability by the pairing available
    panel_distribution(axes[0], rel)
    panel_expression(axes[1], rows)
    panel_pairing(axes[2], rel_by_pairing)

    for ax in axes:
        tidy(ax)

    # the title says why the numbers matter
    fig.suptitle(
        "Every correlation in this analysis is capped by these numbers: "
        "a gene measured once, badly, cannot correlate with anything",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    path = OUT / "ish_reliability.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f"\n{path}")
