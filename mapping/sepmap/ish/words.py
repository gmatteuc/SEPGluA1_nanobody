"""Annotation words and GO terms at the top of the gene ranking, tested.

ish.compare leaves one rho per gene, and the only grouping so far was the
hand-written `category` column of gene_targets.csv. It is coarse ("trafficking"
holds Cacng8 beside the presynaptic Bsn and Syn1), and explaining the ranking with
it is close to circular, since the categories are the reasons the genes were
picked. So the question is asked again with annotation nobody here wrote: each
gene's GO terms from mygene.info, cached one JSON per gene, give two feature sets:

    terms    whole GO terms, such as "regulation of receptor recycling"
    words    single words of those terms and of the gene's name, such as
             "recycling", "endosome", "postsynaptic": coarser, but they pool
             terms that say the same thing differently

For each feature held by at least ish_words.min_genes genes (and by at most
max_share of them), the genes with it are compared with the genes without it on
their rho with the adult nano map: Mann-Whitney on rho, the gap between the two
medians as the effect size with a percentile bootstrap interval,
Benjamini-Hochberg across features. Ranks rather than means, because 95 genes is
small and rho is bounded. Sorted by effect size, a five-gene group comes first for
free; its interval shows it.

Genes sharing a GO term are co-expressed, so their rho values are not independent
draws and the p values are anticonservative by a lot (the null that Fulcher 2021
measured inflated 875-fold): what counts is the effect size and the order, with q
only to separate features of similar size, and the figure says so. The pool is 95
genes chosen for AMPA receptors, so "synaptic words win" partly describes how the
panel was built; contrasts within it, such as postsynaptic against presynaptic
(named before the test was run), are safer.

The bootstrap of each reading draws from one generator, feature after feature in
the order Python sets give them, so the intervals repeat only with PYTHONHASHSEED
fixed.

Run by run_ish_words.py.
"""

import csv
import json
import re
import time
import urllib.parse
import urllib.request
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from scipy.stats import false_discovery_control, mannwhitneyu

from sepmap.config import DATA, SETTINGS
from sepmap.plotting import RED

# the reading of the figure and of the printed summary; the genes a feature needs,
# the largest share it may have, the shortest word, the features shown and the
# bootstrap resamples
ISH = SETTINGS["ish"]
ISH_WORDS = SETTINGS["ish_words"]

RHO = DATA / "adult_v2" / "ish" / "gene_correlations.csv"
OUT = DATA / "adult_v2" / "ish"
CACHE = OUT / "annotation"

MYGENE = "https://mygene.info/v3/query"

# words named before the test was run: if the nanobody reports surface GluA1 at
# glutamatergic postsynapses, postsynaptic words should track the map and
# presynaptic ones should not (the panel holds 68 genes against 49, so its design
# does not force the split)
PREDICTED = (
    "postsynaptic",
    "presynaptic",
    "axon",
    "glutamatergic",
    "dendritic",
    "spine",
    "vesicle",
    "inhibitory",
    "gabaergic",
)

# the contrast drawn in the figure
CONTRAST = ("postsynaptic", "presynaptic")

# words that appear in GO terms as grammar rather than as meaning
STOP = {
    "process",
    "activity",
    "involved",
    "response",
    "positive",
    "negative",
    "protein",
    "cellular",
    "from",
    "into",
    "with",
    "that",
    "this",
    "other",
    "type",
    "like",
    "containing",
    "complex",
}


def fetch(symbol: str) -> dict:
    """The mygene.info record of one mouse gene, cached on disk.

    Cached because it is the only step that needs a server, and a re-run should
    not ask 95 times for an answer that does not change. Delete the folder to
    refresh.
    """
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / (symbol + ".json")
    if path.exists():
        with open(path, encoding="utf-8") as fh:
            return json.load(fh)

    # the record whose symbol matches, or an empty one
    query = urllib.parse.urlencode(
        {
            "q": f"symbol:{symbol}",
            "species": "mouse",
            "fields": "symbol,name,summary,go",
            "size": 5,
        }
    )
    with urllib.request.urlopen(f"{MYGENE}?{query}", timeout=60) as fh:
        hits = json.load(fh).get("hits", [])
    hits = [h for h in hits if h.get("symbol", "").lower() == symbol.lower()]
    record = hits[0] if hits else {}
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(record, fh, indent=1)

    # a courtesy to a free public API
    time.sleep(0.2)
    return record


def go_terms(record: dict) -> set[str]:
    """Every GO term on the record, as a set of strings.

    mygene returns a dict when a gene has one term in a branch and a list when
    it has several, so both shapes are handled. Evidence codes are kept in the
    cached JSON but not filtered on: dropping IEA would remove most of the
    annotation for the less-studied genes and bias the comparison towards the
    ones somebody has studied by hand.
    """
    terms = set()
    for branch in record.get("go", {}).values():
        for item in branch if isinstance(branch, list) else [branch]:
            if isinstance(item, dict) and item.get("term"):
                terms.add(item["term"].strip().lower())
    return terms


def words_of(terms: set[str], name: str | None) -> set[str]:
    """Single words from the GO terms and the gene's full name."""
    out = set()
    for text in list(terms) + [name or ""]:
        for w in re.split(r"[^a-z]+", text.lower()):
            if len(w) < ISH_WORDS["min_word"] or w in STOP:
                continue

            # a very light plural rule: enough to pool "spine" and "spines", held
            # back from "across", "synapsis", "exocytosis" and the like
            if (
                w.endswith("s")
                and len(w) > ISH_WORDS["min_word"] + 1
                and w[-2:] not in ("ss", "us", "is", "as")
            ):
                w = w[:-1]
            out.add(w)
    return out


def load_rho() -> dict[str, dict[str, float]]:
    """{reading: {gene: rho}} from what ish.compare wrote."""
    per = defaultdict(dict)
    with open(RHO, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            per[r["reading"]][r["symbol"]] = float(r["rho"])
    return per


def gap_interval(
    a: np.ndarray, b: np.ndarray, rng: np.random.Generator
) -> tuple[float, float]:
    """A 95% percentile bootstrap interval for the difference of the two medians.

    Without it the panel misleads: a five-gene group reaches a large median gap
    for nothing, and sorted by effect size such groups come first. The interval
    shows that on the figure.
    """
    n_boot = ISH_WORDS["n_boot"]
    da = np.median(rng.choice(a, (n_boot, len(a))), axis=1)
    db = np.median(rng.choice(b, (n_boot, len(b))), axis=1)
    lo, hi = np.percentile(da - db, [2.5, 97.5])
    return float(lo), float(hi)


def test_features(features: dict[str, set[str]], rho: dict[str, float]) -> list[dict]:
    """One row per feature: the two medians, the gap, Mann-Whitney p, BH q.

    `features` is {feature: set of genes}. Genes without a rho for this reading are
    absent from both sides. Rows come sorted by gap, largest first.
    """
    genes = set(rho)
    min_genes = ISH_WORDS["min_genes"]

    # seeded, so a re-run gives the same intervals (with the features in the same
    # order)
    rng = np.random.default_rng(0)
    rows = []
    for feature, carriers in features.items():
        inside = sorted(carriers & genes)
        outside = sorted(genes - carriers)
        if (
            not min_genes <= len(inside) <= ISH_WORDS["max_share"] * len(genes)
            or len(outside) < min_genes
        ):
            continue

        # the genes with the feature against the genes without
        a = np.array([rho[g] for g in inside])
        b = np.array([rho[g] for g in outside])
        _, p = mannwhitneyu(a, b, alternative="two-sided")
        lo, hi = gap_interval(a, b, rng)
        rows.append(
            dict(
                feature=feature,
                n_genes=len(inside),
                median_with=float(np.median(a)),
                median_without=float(np.median(b)),
                gap=float(np.median(a) - np.median(b)),
                gap_lo=lo,
                gap_hi=hi,
                p=float(p),
                genes=" ".join(inside),
            )
        )

    # Benjamini-Hochberg across the features
    if rows:
        q = false_discovery_control([r["p"] for r in rows])
        for r, qi in zip(rows, q):
            r["q"] = float(qi)
    return sorted(rows, key=lambda r: -r["gap"])


def bars(ax: Axes, rows: list[dict], title: str, xlim: tuple[float, float]) -> None:
    """Draw the top features by effect size as bars; the lower q, the darker."""
    sel = rows[: ISH_WORDS["n_shown"]][::-1]
    y = np.arange(len(sel))
    shade = [str(np.clip(0.75 * min(r["q"], 1.0) ** 0.5, 0.05, 0.8)) for r in sel]
    gaps = np.array([r["gap"] for r in sel])
    err = np.abs(
        np.array(
            [[r["gap"] - r["gap_lo"] for r in sel], [r["gap_hi"] - r["gap"] for r in sel]]
        )
    )
    ax.barh(y, gaps, color=shade, edgecolor="0.25", linewidth=0.5)
    ax.errorbar(gaps, y, xerr=err, fmt="none", ecolor="0.2", elinewidth=0.8, capsize=2)
    ax.set_yticks(y)
    ax.set_yticklabels(
        [f"{r['feature'][:46]}  ({r['n_genes']})" for r in sel], fontsize=7
    )
    ax.axvline(0, color="0.3", lw=0.6)
    ax.set_xlim(xlim)
    ax.set_xlabel("median rho with the feature  -  median rho without", fontsize=8)
    ax.set_title(title, fontsize=9)
    ax.tick_params(axis="x", labelsize=7)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def strip(ax: Axes, rows: list[dict], rho: dict[str, float]) -> None:
    """Draw the contrast named in advance as the distributions behind it.

    Not whichever feature came first: that one is usually a five-gene group, and
    picking it after the fact is how a ranked list turns into a claim.
    """
    by = {r["feature"]: r for r in rows}
    rng = np.random.default_rng(0)
    labels, i = [], 0
    for word in CONTRAST:
        r = by.get(word)
        if r is None:
            continue
        carriers = set(r["genes"].split())
        for has, colour in ((False, "0.7"), (True, RED)):
            values = [rho[g] for g in sorted(carriers if has else set(rho) - carriers)]
            ax.scatter(
                np.full(len(values), i) + rng.uniform(-0.11, 0.11, len(values)),
                values,
                s=13,
                facecolor=colour,
                edgecolor="0.25",
                linewidth=0.4,
                zorder=2,
            )
            ax.plot(
                [i - 0.28, i + 0.28],
                [np.median(values)] * 2,
                color="0.15",
                lw=1.6,
                zorder=3,
            )
            labels.append(f"{'yes' if has else 'no'}\n{len(values)}")
            i += 1
        ax.text(
            i - 1.5,
            ax.get_ylim()[0],
            f"{word}\ngap {r['gap']:+.3f}, q = {r['q']:.3f}",
            ha="center",
            va="bottom",
            fontsize=7.5,
        )
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=7.5)
    ax.set_ylabel(f"rho with the adult nano map ({ISH['reading']})", fontsize=8)
    ax.set_title("the contrast predicted in advance", fontsize=9)
    ax.tick_params(axis="y", labelsize=7)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def figure(terms: list[dict], words: list[dict], rho: dict[str, float]) -> None:
    """Draw the two bar panels and the contrast; saved as ish_word_enrichment.png."""
    fig, axes = plt.subplots(
        1, 3, figsize=(15.5, 4.6), gridspec_kw=dict(width_ratios=[1.3, 1.0, 0.8])
    )

    # one x range for both bar panels, from the intervals of the features shown
    n_shown = ISH_WORDS["n_shown"]
    xlim = (
        min(-0.02, min(r["gap_lo"] for r in terms[:n_shown] + words[:n_shown]) - 0.02),
        max(r["gap_hi"] for r in terms[:n_shown] + words[:n_shown]) + 0.02,
    )
    bars(axes[0], terms, "GO terms", xlim)
    bars(axes[1], words, "words, from GO terms and gene names", xlim)
    strip(axes[2], words, rho)

    # the title says how to read the bars
    fig.suptitle(
        "What the top of the gene ranking is made of.  Bars are ordered by effect "
        "size, whiskers are a 95% bootstrap interval, darkness is the BH q."
        "\nThe p behind that q is anticonservative, because genes sharing a term are "
        "co-expressed (Fulcher 2021): read the size and the width, not the significance.",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    path = OUT / "ish_word_enrichment.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f"\n{path}")


def feature_genes(genes: list[str]) -> tuple[dict[str, set[str]], dict[str, set[str]]]:
    """The genes of each GO term and of each word; prints how many of each."""
    term_genes, word_genes, no_go = defaultdict(set), defaultdict(set), []
    for gene in genes:
        record = fetch(gene)
        terms = go_terms(record)
        if not terms:
            no_go.append(gene)
        for t in terms:
            term_genes[t].add(gene)
        for w in words_of(terms, record.get("name", "")):
            word_genes[w].add(gene)
    print(
        f"{len(term_genes)} distinct GO terms, {len(word_genes)} distinct words"
        + (f"; no GO annotation for {no_go}" if no_go else "")
    )
    return term_genes, word_genes


def enrichment_rows(
    rho: dict[str, dict[str, float]],
    term_genes: dict[str, set[str]],
    word_genes: dict[str, set[str]],
) -> list[dict]:
    """Every feature tested, for every reading: the rows of feature_enrichment.csv."""
    rows = []
    for reading in sorted(rho):
        for kind, features in (("term", term_genes), ("word", word_genes)):
            for r in test_features(features, rho[reading]):
                rows.append(dict(reading=reading, kind=kind, **r))
    return rows


def write_table(rows: list[dict]) -> None:
    """Write feature_enrichment.csv, floats to four significant digits."""
    path = OUT / "feature_enrichment.csv"
    fields = [
        "reading",
        "kind",
        "feature",
        "n_genes",
        "median_with",
        "median_without",
        "gap",
        "gap_lo",
        "gap_hi",
        "p",
        "q",
        "genes",
    ]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(
            {k: (f"{v:.4g}" if isinstance(v, float) else v) for k, v in r.items()}
            for r in rows
        )
    print(f"{len(rows)} rows -> {path}")


def report(rows: list[dict], plot_reading: str) -> tuple[list[dict], list[dict]]:
    """Print the words named in advance, then the top features, for the reading shown.

    Returns that reading's rows of GO terms and of words.
    """
    terms = [r for r in rows if r["reading"] == plot_reading and r["kind"] == "term"]
    words = [r for r in rows if r["reading"] == plot_reading and r["kind"] == "word"]
    by = {r["feature"]: r for r in words}
    print(f"\nthe words named in advance, for {plot_reading}")
    for word in PREDICTED:
        r = by.get(word)
        if r:
            result = f"gap {r['gap']:+.3f}  n={r['n_genes']:3d}  q={r['q']:.3f}"
        else:
            result = "not tested (too few genes, or absent)"
        print(f"  {word:15s} " + result)
    for label, sel in (("GO terms", terms), ("words", words)):
        print(f"\ntop {label} for {plot_reading}  (gap in median rho, n genes, BH q)")
        for r in sel[: ISH_WORDS["n_shown"]]:
            print(
                f"  {r['gap']:+.3f} [{r['gap_lo']:+.2f} {r['gap_hi']:+.2f}]  "
                f"n={r['n_genes']:3d}  q={r['q']:.3f}  {r['feature']}"
            )
    return terms, words


def main() -> None:
    """Test every GO term and word against the gene ranking, write, report and draw."""
    # the rho of each gene and reading
    plot_reading = ISH["reading"]
    rho = load_rho()
    genes = sorted(rho[plot_reading])
    print(f"{len(genes)} genes, readings {sorted(rho)}")

    # the genes of each GO term and of each word, every feature tested
    term_genes, word_genes = feature_genes(genes)
    rows = enrichment_rows(rho, term_genes, word_genes)
    write_table(rows)

    # the words named in advance, then the top features, for the reading shown
    terms, words = report(rows, plot_reading)
    figure(terms, words, rho[plot_reading])
