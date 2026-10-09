"""One gene table for every ISH analysis: both panels, their experiments, QC and labels.

Two panels were built at different times for different questions: P9's 100
hand-picked genes (mapping/gene_targets.csv, one coronal experiment each) and the
390-gene ontology panel (ish.panel_build, every experiment of genes taken from GO
terms). Each analysis read one or the other. Here they become one table, a row per
experiment:

    experiments     the union of the two panels (numbers_gene_table.csv counts
                    them). A P9 gene whose own grid cannot be used (Allen returns
                    404 for Chat and Tph2; Olig2 and Calb2 sit in a box of their
                    own) and that the ontology panel lacks is repaired with its
                    other Allen experiments, fetched once (the repair). Sst needs
                    none: the ontology panel holds two of its experiments
    region means    per experiment and structure, full and eroded by one 200 um
                    voxel (ish.regions), with the sections ish.section_qc flagged
                    set missing first
    reliability     per gene with two experiments or more, the median Spearman
                    between them (ish.reliability)
    profile         per gene, the mean of its experiments' ranks (reliability.merge),
                    the one profile every analysis correlates with the map
    labels          the ontology panel's role, P9's category (a label only, never a
                    group of an analysis), the gene's full name and GO annotation
                    from mygene.info, and the gene sets of analysis 3
                    (ish.gene_sets), with GO's ancestry from a dated go-basic.obo

An experiment that cannot be used keeps its row, excluded, with the reason. Network
answers are cached under adult_v2/ish_analysis/cache/: the Allen experiment lists of
the repaired genes, the mygene records (those of 5 October in adult_v2/ish/annotation/
are read first and never rewritten) and go-basic.obo with its release.

gene_documentation.csv is the table for the supplement Sami asked for (28 April):
per gene, what it is, why it is in the panel, which experiments were used, how
reliable its map is; a reference column is left for curation by hand. It is written
UTF-8 with a byte-order mark so that Excel opens it as it is.

Run by run_ish_section_qc.py (the experiment list) and run_ish_gene_table.py.
"""

import json
import re
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from sepmap.config import DATA, SETTINGS
from sepmap.ish import panel_build, panel_fetch, regions, reliability, section_qc
from sepmap.ish.gene_sets import context_set, set_members
from sepmap.ish.numbers import numbers_frame, numbers_path
from sepmap.structures import ISH_OUT, TABLES
from sepmap.volumes.per_mouse import structure_terms

# the voxels and structures a gene value needs to be used; the reliability below
# which a gene's Allen map does not reproduce
ISH = SETTINGS["ish"]
ISH_PANEL_TEST = SETTINGS["ish_panel_test"]

CACHE = ISH_OUT / "cache"
EXPERIMENTS = TABLES / "experiments.csv"
REGION_TABLE = TABLES / "gene_region_table.csv"
GENE_TABLE = TABLES / "gene_table.csv"
PROFILES = TABLES / "gene_profiles.csv"
DOCUMENTATION = TABLES / "gene_documentation.csv"
NUMBERS = numbers_path("gene_table")

# the mygene records of 5 October, read first and never written
OLD_MYGENE = DATA / "adult_v2" / "ish" / "annotation"

MYGENE = "https://mygene.info/v3/query"
GO_OBO = "http://purl.obolibrary.org/obo/go/go-basic.obo"

# genes asked of mygene.info per request
MYGENE_BATCH = 200

# how the downloads introduce themselves
USER_AGENT = "sepmap (SEP-GluA1 histology analysis)"

# the genes whose reliability the text quotes: the gene of the stained protein, the
# gene named for the leftover, and a scaffold near the top of the ranking
QUOTED_GENES = ("Gria1", "Cacng8", "Dlg2")


class OfflineError(RuntimeError):
    """A cache lacks an answer and the run was told not to ask the server."""


# ===== The experiments =====


def read_panel(name: str) -> pd.DataFrame:
    """One panel of settings.toml [ish_panels], experiment ids as strings."""
    path, _ = regions.panel_files(name)
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def grid_problem(experiment_id: str) -> str:
    """Why an experiment's grid cannot be used, or '' when it can."""
    mhd = regions.ISH_DIR / f"{experiment_id}_energy.mhd"
    if not panel_fetch.already_there(experiment_id):
        return "grid not on disk"
    dims = panel_fetch.dims_of(mhd)
    if dims != regions.GRID_DIMS:
        return f"grid is {dims}, not the reference box {regions.GRID_DIMS}"
    return ""


def repair_genes(p9: pd.DataFrame, ontology: pd.DataFrame) -> list[str]:
    """P9 genes whose own grid cannot be used and that the ontology panel lacks."""
    in_ontology = set(ontology["symbol"])
    out = []
    for _, r in p9.iterrows():
        if r["symbol"] not in in_ontology and grid_problem(r["experiment_id"]):
            out.append(r["symbol"])
    return out


def repair_rows(genes: list[str], p9: pd.DataFrame, offline: bool) -> list[dict]:
    """The other Allen experiments of the genes to repair, their grids fetched.

    The experiment lists are cached under cache/allen/; a grid already on disk is
    left alone, and one that cannot be fetched keeps its row, so the gene table
    says why it is not used.
    """
    allen_cache = CACHE / "allen"
    own = set(p9["experiment_id"])
    rows = []
    for symbol in genes:
        if offline and not (allen_cache / f"allen_{symbol}.json").exists():
            raise OfflineError(
                f"no cached Allen experiment list for {symbol} in {allen_cache}; "
                "run once without --offline"
            )
        for e in panel_build.allen_experiments(symbol, cache=allen_cache):
            eid = str(e["id"])
            if eid in own:
                continue
            if not panel_fetch.already_there(eid) and not offline:
                why = panel_fetch.fetch(eid)
                if why:
                    print(f"  repair: {symbol} {eid} not fetched ({why})", flush=True)
            rows.append(dict(symbol=symbol, experiment_id=eid, plane=e["plane"]))
    return rows


def union_rows(
    p9: pd.DataFrame, ontology: pd.DataFrame, repair: list[dict]
) -> pd.DataFrame:
    """One row per experiment of the two panels and of the repair, with its sources.

    An experiment both panels list is one row. Columns: symbol, experiment_id,
    plane, p9_experiment, ontology_experiment, repair_experiment, sorted by gene
    and experiment.
    """
    rows = {}
    for panel in (p9, ontology):
        for _, r in panel.iterrows():
            rows.setdefault(
                r["experiment_id"],
                dict(
                    symbol=r["symbol"], experiment_id=r["experiment_id"], plane=r["plane"]
                ),
            )
    for r in repair:
        rows.setdefault(r["experiment_id"], r)
    table = pd.DataFrame(list(rows.values()))
    table["p9_experiment"] = table["experiment_id"].isin(set(p9["experiment_id"]))
    table["ontology_experiment"] = table["experiment_id"].isin(
        set(ontology["experiment_id"])
    )
    table["repair_experiment"] = ~table["p9_experiment"] & ~table["ontology_experiment"]
    return table.sort_values(["symbol", "experiment_id"]).reset_index(drop=True)


def experiment_list(offline: bool = False) -> pd.DataFrame:
    """Every experiment of the two panels and of the repair (union_rows)."""
    p9 = read_panel("targets")
    ontology = read_panel("ontology")
    repaired = repair_genes(p9, ontology)
    return union_rows(p9, ontology, repair_rows(repaired, p9, offline))


def load_experiments() -> pd.DataFrame:
    """The experiment list run_ish_section_qc wrote."""
    if not EXPERIMENTS.exists():
        raise FileNotFoundError(
            f"{EXPERIMENTS} not found: run run_ish_section_qc.py first"
        )
    table = pd.read_csv(EXPERIMENTS, dtype={"experiment_id": str})
    for column in ("p9_experiment", "ontology_experiment", "repair_experiment"):
        table[column] = table[column].astype(str) == "True"
    return table


# ===== Region means, reliability, profiles =====


def region_rows(
    experiments: pd.DataFrame,
    flags: dict[str, list[int]],
    ann: np.ndarray,
    names: dict[int, str],
    eroded: np.ndarray,
) -> list[dict]:
    """Full and eroded mean per experiment and structure, flagged sections missing.

    `flags` holds the sections to set missing per experiment id (empty: no QC).
    Experiments whose grid cannot be used are skipped; the gene table says why.
    """
    rows = []
    for i, r in enumerate(experiments.to_dict("records"), 1):
        eid = r["experiment_id"]
        if grid_problem(eid):
            continue
        vol = section_qc.read_grid(eid, ann.shape)
        axis = section_qc.SECTION_AXIS[r["plane"]]
        vol = section_qc.apply_flags(vol, flags.get(eid, []), axis)
        means = regions.region_means(vol, ann, names, eroded)
        for structure, (full, ero, n, n_total) in means.items():
            rows.append(
                dict(
                    symbol=r["symbol"],
                    experiment_id=eid,
                    plane=r["plane"],
                    structure=structure,
                    ish_mean=full,
                    ish_mean_eroded=ero,
                    n_voxels=n,
                    n_voxels_structure=n_total,
                    coverage=n / max(n_total, 1),
                )
            )
        if i % 100 == 0:
            print(f"  {i}/{len(experiments)} experiments", flush=True)
    return rows


def region_table(experiments: pd.DataFrame, flags: dict[str, list[int]]) -> pd.DataFrame:
    """region_rows on the CCF of the 200 um grid, full and eroded, as one table."""
    ann = regions.annotation_200()
    names, _, _ = structure_terms()
    eroded = regions.eroded_annotation(ann)
    return pd.DataFrame(region_rows(experiments, flags, ann, names, eroded))


def experiment_profiles(
    region: pd.DataFrame, column: str = "ish_mean"
) -> dict[str, dict[str, dict[str, float]]]:
    """{gene: {experiment: {structure: energy}}}, with ish.min_voxels voxels of data.

    `column` is the mean taken, ish_mean or ish_mean_eroded; a structure whose
    eroded mean is missing (too small to erode) is left out of that profile.
    """
    usable = region[(region["n_voxels"] >= ISH["min_voxels"]) & region[column].notna()]
    per = defaultdict(lambda: defaultdict(dict))
    for symbol, eid, structure, value in zip(
        usable["symbol"], usable["experiment_id"], usable["structure"], usable[column]
    ):
        per[symbol][eid][structure] = float(value)
    return per


def gene_reliability(per: dict) -> pd.DataFrame:
    """Per gene: experiments, pairs, the median Spearman between them, the level.

    The level is the median energy over every structure of every experiment, as
    ish.reliability takes it.
    """
    rows = []
    for symbol, profiles in sorted(per.items()):
        pairs = reliability.pair_reliability(profiles)
        rel = float(np.median([p[2] for p in pairs])) if pairs else np.nan
        level = float(np.median([v for p in profiles.values() for v in p.values()]))
        rows.append(
            dict(
                symbol=symbol,
                n_experiments_used=len(profiles),
                n_pairs=len(pairs),
                reliability=rel,
                median_energy=level,
            )
        )
    return pd.DataFrame(rows)


def used_experiment_profiles(
    table: pd.DataFrame,
) -> dict[str, dict[str, dict[str, float]]]:
    """experiment_profiles of the experiments the gene table uses, from the region table.

    `table` is the gene table; {gene: {experiment: {structure: energy}}}.
    """
    region = load_region_table()
    used = set(table.loc[~table["excluded"], "experiment_id"])
    return experiment_profiles(region[region["experiment_id"].isin(used)])


def gene_levels(table: pd.DataFrame) -> tuple[dict[str, float], dict[str, float]]:
    """{gene: reliability} and {gene: median energy}, from the usable experiments.

    `table` is the gene table; a gene measured once has no reliability.
    """
    rel = gene_reliability(used_experiment_profiles(table)).set_index("symbol")
    return rel["reliability"].dropna().to_dict(), rel["median_energy"].to_dict()


def profile_rows(per: dict) -> list[dict]:
    """Each gene's merged profile: the mean of its experiments' 0-1 ranks."""
    rows = []
    for symbol, profiles in sorted(per.items()):
        merged, n_exp = reliability.merge(profiles)
        for structure in sorted(merged):
            rows.append(
                dict(
                    symbol=symbol,
                    structure=structure,
                    rank_mean=merged[structure],
                    n_experiments=n_exp[structure],
                )
            )
    return rows


def load_region_table() -> pd.DataFrame:
    """The region means run_ish_gene_table wrote, flagged sections already missing."""
    if not REGION_TABLE.exists():
        raise FileNotFoundError(
            f"{REGION_TABLE} not found: run run_ish_gene_table.py first"
        )
    return pd.read_csv(REGION_TABLE, dtype={"experiment_id": str})


def load_profiles() -> dict[str, dict[str, float]]:
    """{gene: {structure: mean rank}}, from gene_profiles.csv."""
    if not PROFILES.exists():
        raise FileNotFoundError(f"{PROFILES} not found: run run_ish_gene_table.py first")
    table = pd.read_csv(PROFILES)
    out = defaultdict(dict)
    for symbol, structure, value in zip(
        table["symbol"], table["structure"], table["rank_mean"]
    ):
        out[symbol][structure] = float(value)
    return dict(out)


# ===== Annotation =====


def cached_record(symbol: str) -> dict | None:
    """The mygene record of a gene from a cache, or None when neither has it."""
    for folder in (OLD_MYGENE, CACHE / "mygene"):
        path = folder / f"{symbol}.json"
        if path.exists():
            with open(path, encoding="utf-8") as fh:
                return json.load(fh)
    return None


def fetch_records(symbols: list[str]) -> None:
    """Ask mygene.info for the records of `symbols`, MYGENE_BATCH at a time; cache them.

    The record kept is the first hit whose symbol matches exactly, as words.fetch
    keeps it; a gene with none is cached as an empty record, so it is not asked
    again.
    """
    folder = CACHE / "mygene"
    folder.mkdir(parents=True, exist_ok=True)
    for start in range(0, len(symbols), MYGENE_BATCH):
        batch = symbols[start : start + MYGENE_BATCH]
        data = urllib.parse.urlencode(
            {
                "q": ",".join(batch),
                "scopes": "symbol",
                "species": "mouse",
                "fields": "symbol,name,summary,go",
            }
        ).encode()
        request = urllib.request.Request(MYGENE, data=data)
        with urllib.request.urlopen(request, timeout=120) as fh:
            hits = json.load(fh)
        for symbol in batch:
            match = [
                h
                for h in hits
                if h.get("query") == symbol
                and h.get("symbol", "").lower() == symbol.lower()
            ]
            record = {k: v for k, v in match[0].items() if k != "query"} if match else {}
            with open(folder / f"{symbol}.json", "w", encoding="utf-8") as fh:
                json.dump(record, fh, indent=1)

        # a courtesy to a free public API
        time.sleep(0.5)


def mygene_records(symbols: list[str], offline: bool) -> dict[str, dict]:
    """The mygene record of every gene, from the caches or, for the rest, the server."""
    missing = [s for s in symbols if cached_record(s) is None]
    if missing and offline:
        raise OfflineError(
            f"{len(missing)} genes have no cached mygene record ({', '.join(missing[:5])}"
            "...); run once without --offline"
        )
    if missing:
        print(f"  asking mygene.info for {len(missing)} genes", flush=True)
        fetch_records(missing)
    return {s: cached_record(s) or {} for s in symbols}


def go_annotations(record: dict, branch: str) -> set[str]:
    """The GO ids of one branch (BP, CC or MF) of a record, NOT qualifiers left out.

    mygene gives a dict when a gene has one term in a branch and a list when it
    has several.
    """
    items = record.get("go", {}).get(branch, [])
    if isinstance(items, dict):
        items = [items]
    out = set()
    for item in items:
        if "NOT" in item.get("qualifier", ""):
            continue
        if item.get("id"):
            out.add(item["id"])
    return out


def load_obo(offline: bool) -> tuple[dict[str, set[str]], dict[str, str], str]:
    """GO's parents by is_a and part_of, term names, and the release, from go-basic.obo.

    Downloaded once into the cache; go-basic holds only relations that never cross
    the three branches, so ancestry stays within cellular component.
    """
    path = CACHE / "go-basic.obo"
    if not path.exists():
        if offline:
            raise OfflineError(f"{path} not cached; run once without --offline")
        CACHE.mkdir(parents=True, exist_ok=True)
        print("  downloading go-basic.obo", flush=True)
        # the Gene Ontology's server refuses Python's default user agent (HTTP 403)
        request = urllib.request.Request(GO_OBO, headers={"User-Agent": USER_AGENT})
        with urllib.request.urlopen(request, timeout=300) as fh:
            path.write_bytes(fh.read())
    return read_obo(path)


def read_obo(path: Path) -> tuple[dict[str, set[str]], dict[str, str], str]:
    """GO's parents by is_a and part_of, term names, and the release, from an obo file.

    An alt id's parent is its primary term, so an annotation to an alt id has the
    primary term's ancestry.
    """
    parents, names, release = defaultdict(set), {}, ""
    term = None
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("data-version:") and not release:
                release = line.split(":", 1)[1].strip()
            if line.startswith("["):
                term = None
            if line == "[Term]":
                term = ""
            elif term is not None and line.startswith("id: "):
                term = line[4:]
            elif term and line.startswith("name: "):
                names[term] = line[6:]
            elif term and line.startswith("alt_id: "):
                parents[line[8:]].add(term)
            elif term and line.startswith("is_a: "):
                parents[term].add(line[6:].split(" ! ")[0].strip())
            elif term and line.startswith("relationship: part_of "):
                parents[term].add(line.split()[2])
    return dict(parents), names, release


def with_go_names(text: str, names: dict[str, str]) -> str:
    """`text` with each GO id followed by its name in brackets, unless named already."""

    def named(match: re.Match) -> str:
        """The GO id matched, with its name after it when the text lacks it."""
        term = match.group(0)
        after = text[match.end() : match.end() + 1 + len(names.get(term, ""))]
        if term not in names or after.strip() == names[term]:
            return term
        return f"{term} ({names[term]})"

    return re.sub(r"GO:\d{7}", named, text)


def go_term_list(terms: list[str], names: dict[str, str]) -> str:
    """GO ids with their names, joined: 'GO:0014069 postsynaptic density; ...'."""
    return "; ".join(f"{t} {names.get(t, '')}".strip() for t in terms)


def with_ancestors(terms: set[str], parents: dict[str, set[str]]) -> set[str]:
    """The terms and every ancestor of them, by the relations in `parents`."""
    out = set()
    todo = list(terms)
    while todo:
        t = todo.pop()
        if t in out:
            continue
        out.add(t)
        todo.extend(parents.get(t, ()))
    return out


# ===== The tables =====


def cellular_components(
    symbols: list[str], records: dict[str, dict], parents: dict[str, set[str]]
) -> dict[str, set[str]]:
    """{gene: its GO cellular-component terms with their ancestors}, from its record."""
    return {
        g: with_ancestors(go_annotations(records.get(g, {}), "CC"), parents)
        for g in symbols
    }


def gene_labels(
    experiments: pd.DataFrame, records: dict[str, dict], parents: dict[str, set[str]]
) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    """Per gene: P9's category and description, the panel role, name, sets.

    Returns the per-gene table and the members of each set (ish.gene_sets) among
    the genes with at least one usable experiment.
    """
    p9 = read_panel("targets").drop_duplicates("symbol").set_index("symbol")
    ontology = read_panel("ontology").drop_duplicates("symbol").set_index("symbol")
    genes = sorted(experiments["symbol"].unique())
    role = {g: ontology["role"].get(g, "") for g in genes}
    components = cellular_components(genes, records, parents)
    usable = sorted(experiments.loc[~experiments["excluded"], "symbol"].unique())
    members = set_members(usable, role, components)
    in_sets = defaultdict(list)
    for name, symbols in members.items():
        for s in symbols:
            in_sets[s].append(name)
    rows = []
    for g in genes:
        record = records.get(g, {})
        rows.append(
            dict(
                symbol=g,
                name=record.get("name", ""),
                p9_gene=g in p9.index,
                p9_category=p9["category"].get(g, ""),
                p9_description=p9["description"].get(g, ""),
                ontology_role=role[g],
                ontology_go_terms=ontology["go_terms"].get(g, ""),
                gene_sets="; ".join(in_sets.get(g, [])),
                n_go_cc=len(go_annotations(record, "CC")),
                summary=record.get("summary", ""),
            )
        )
    return pd.DataFrame(rows), members


def context_members(symbols: list[str], offline: bool = True) -> list[str]:
    """The genes of the context group of ish.gene_sets, from the cached GO records.

    The same records, ontology and roles as gene_labels; with `offline` a gene
    missing from the caches stops the run rather than asking the server.
    """
    ontology = read_panel("ontology").drop_duplicates("symbol").set_index("symbol")
    role = {g: ontology["role"].get(g, "") for g in symbols}
    records = mygene_records(symbols, offline)
    parents, _, _ = load_obo(offline)
    return context_set(symbols, role, cellular_components(symbols, records, parents))


def exclusion_reasons(experiments: pd.DataFrame) -> pd.DataFrame:
    """The experiment list with excluded and exclude_reason, for grids not usable.

    A grid not on disk carries the reason ish.panel_fetch recorded, when it did.
    """
    failed = regions.fetch_failures()
    out = experiments.copy()
    reasons = []
    for eid in out["experiment_id"]:
        why = grid_problem(eid)
        if why == "grid not on disk" and eid in failed:
            why += f" ({failed[eid]})"
        reasons.append(why)
    out["exclude_reason"] = reasons
    out["excluded"] = out["exclude_reason"] != ""
    return out


def coverage_exclusions(experiments: pd.DataFrame, region: pd.DataFrame) -> pd.DataFrame:
    """The experiment list with the grids that measure too few structures excluded.

    An experiment enters a gene's profile only with ish.min_structures structures of
    ish.min_voxels voxels with data, the number a correlation needs: the merged
    profile averages each experiment's ranks put on 0-1, so one that measures two
    structures would give them rank 0 and 1 and pull the gene's profile there. A
    few early sagittal grids, and three with no data inside the brain, fall here.
    """
    usable = region[region["n_voxels"] >= ISH["min_voxels"]]
    n = usable.groupby("experiment_id")["structure"].nunique()
    out = experiments.copy()
    for i, r in out.iterrows():
        if r["excluded"]:
            continue
        n_structures = int(n.get(r["experiment_id"], 0))
        if n_structures == 0:
            why = "no data inside the brain"
        elif n_structures < ISH["min_structures"]:
            why = (
                f"data in {n_structures} structures, fewer than the "
                f"{ISH['min_structures']} a correlation needs"
            )
        else:
            continue
        out.loc[i, "excluded"] = True
        out.loc[i, "exclude_reason"] = why
    return out


def experiment_table(
    experiments: pd.DataFrame,
    qc_summary: pd.DataFrame,
    region: pd.DataFrame,
    labels: pd.DataFrame,
    rel: pd.DataFrame,
) -> pd.DataFrame:
    """gene_table.csv: one row per experiment, its QC and use, and its gene's labels."""
    qc = qc_summary.set_index("experiment_id")
    usable = region[region["n_voxels"] >= ISH["min_voxels"]]
    n_structures = usable.groupby("experiment_id")["structure"].nunique()
    level = usable.groupby("experiment_id")["ish_mean"].median()
    out = experiments.copy()
    for column in (
        "n_judged",
        "n_flagged",
        "flagged_sections",
        "n_absence_kept",
        "absence_sections",
    ):
        out[column] = out["experiment_id"].map(qc[column])
    out["n_structures"] = out["experiment_id"].map(n_structures).fillna(0).astype(int)
    out["median_energy"] = out["experiment_id"].map(level)
    gene = labels.set_index("symbol")
    for column in ("p9_gene", "p9_category", "ontology_role", "gene_sets"):
        out[column] = out["symbol"].map(gene[column])
    out["gene_reliability"] = out["symbol"].map(rel.set_index("symbol")["reliability"])
    return out


def load_gene_table() -> pd.DataFrame:
    """The experiment rows run_ish_gene_table wrote, its flags back to booleans."""
    if not GENE_TABLE.exists():
        raise FileNotFoundError(
            f"{GENE_TABLE} not found: run run_ish_gene_table.py first"
        )
    table = pd.read_csv(
        GENE_TABLE,
        dtype={"experiment_id": str},
        keep_default_na=False,
        na_values=[""],
    )
    for column in (
        "p9_experiment",
        "ontology_experiment",
        "repair_experiment",
        "excluded",
        "p9_gene",
    ):
        table[column] = table[column].astype(str) == "True"
    for column in ("p9_category", "ontology_role", "gene_sets", "flagged_sections"):
        table[column] = table[column].fillna("")
    return table


def per_gene(table: pd.DataFrame) -> pd.DataFrame:
    """One row per gene of the gene table: its labels and the experiments it uses.

    Columns: symbol, p9_gene, p9_category, ontology_role, gene_sets, reliability,
    n_experiments_used, experiments_used (ids joined by spaces) and
    p9_experiment_id (P9's own experiment when it is used, else '').
    """
    rows = []
    for symbol, mine in table.groupby("symbol", sort=True):
        used = mine[~mine["excluded"]]
        own = used.loc[used["p9_experiment"], "experiment_id"]
        first = mine.iloc[0]
        rows.append(
            dict(
                symbol=symbol,
                p9_gene=bool(first["p9_gene"]),
                p9_category=first["p9_category"],
                ontology_role=first["ontology_role"],
                gene_sets=first["gene_sets"],
                reliability=first["gene_reliability"],
                n_experiments_used=len(used),
                experiments_used=" ".join(used["experiment_id"]),
                p9_experiment_id=own.iloc[0] if len(own) else "",
            )
        )
    return pd.DataFrame(rows)


def documentation_table(
    labels: pd.DataFrame, experiments: pd.DataFrame, rel: pd.DataFrame
) -> pd.DataFrame:
    """gene_documentation.csv: one row per gene, for the supplement."""
    rel = rel.set_index("symbol")
    rows = []
    for _, g in labels.iterrows():
        mine = experiments[experiments["symbol"] == g["symbol"]]
        used = mine[~mine["excluded"]]
        dropped = mine[mine["excluded"]]
        rows.append(
            dict(
                symbol=g["symbol"],
                name=g["name"],
                why_included=g["p9_description"]
                or f"ontology panel: {g['ontology_role']}",
                p9_category=g["p9_category"],
                ontology_role=g["ontology_role"],
                ontology_go_terms=g["ontology_go_terms"],
                gene_sets=g["gene_sets"],
                experiments_used=" ".join(
                    f"{e} ({p})" for e, p in zip(used["experiment_id"], used["plane"])
                ),
                experiments_not_used=" ".join(
                    f"{e} ({r})"
                    for e, r in zip(dropped["experiment_id"], dropped["exclude_reason"])
                ),
                sections_set_missing=" ".join(
                    f"{e}: {s}"
                    for e, s in zip(used["experiment_id"], used["flagged_sections"])
                    if s
                ),
                reliability=rel["reliability"].get(g["symbol"], np.nan),
                summary=g["summary"],
                reference="",
            )
        )
    return pd.DataFrame(rows)


def numbers_table(
    experiments: pd.DataFrame,
    rel: pd.DataFrame,
    members: dict[str, list[str]],
    labels: pd.DataFrame,
    release: str,
) -> pd.DataFrame:
    """numbers_gene_table.csv: the numbers of this step that the text quotes."""
    used = experiments[~experiments["excluded"]]
    values = rel["reliability"].dropna()
    by_gene = rel.set_index("symbol")["reliability"]
    p9_genes = set(labels.loc[labels["p9_gene"], "symbol"])
    low = ISH_PANEL_TEST["min_reliability"]
    rows = [
        ("genes_listed", experiments["symbol"].nunique(), "genes of both panels"),
        ("experiments_listed", len(experiments), "experiments, repair included"),
        ("experiments_used", len(used), "experiments with a usable grid"),
        ("genes_with_profile", rel["symbol"].nunique(), "genes with a profile"),
        (
            "p9_genes_with_profile",
            len(p9_genes & set(rel["symbol"])),
            f"of P9's {len(p9_genes)} genes",
        ),
        ("genes_reliability", len(values), "genes measured more than once"),
        ("reliability_median", round(float(values.median()), 3), "median reliability"),
        ("reliability_q1", round(float(values.quantile(0.25)), 3), "its first quartile"),
        ("reliability_q3", round(float(values.quantile(0.75)), 3), "its third quartile"),
        (f"reliability_below_{low}", int((values < low).sum()), f"genes below {low}"),
        ("go_release", release, "go-basic.obo release of the ancestry"),
    ]
    for gene in QUOTED_GENES:
        value = by_gene.get(gene, np.nan)
        rows.append(
            (f"reliability_{gene}", round(float(value), 3), f"{gene}'s reliability")
        )
    for name, symbols in members.items():
        rows.append((f"set_{name}", len(symbols), f"genes in the set {name}"))
    for reason, n in (
        experiments.loc[experiments["excluded"], "exclude_reason"].value_counts().items()
    ):
        rows.append(("excluded", int(n), reason))
    return numbers_frame(rows)
