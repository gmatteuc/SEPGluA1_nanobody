"""The synapse-density term of analysis 4: postsynaptic genes chosen blind to the map.

Analysis 4 asks how much of the nano map Gria1 and synapse density leave, so its
density term must not be chosen by how well it predicts the map. It is chosen here
by a rule that reads no nano value, and the rule and the genes it gave were
committed before the model they feed was run (git log shows the order):

    pool        genes of the gene table with density_markers.min_experiments usable
                Allen experiments or more (after section QC), measured
                (ish.min_voxels) in density_markers.min_coverage of the declared
                structures where Gria1 is measured, and annotated in mouse GO to
                the postsynaptic density or the postsynaptic specialization, or to a
                term below either by is_a or part_of (POOL_TERMS)
    excluded    the genes that place or regulate AMPA receptors at the surface
                belong to the surface side, which the leftover may hold, not to
                density: Gria1 to Gria4, the localisation set and the AMPA receptor
                complex family (ish.gene_sets), and every gene annotated to a term of
                EXCLUSION_TERMS or below it, whatever the evidence. Broad plasticity
                terms do not exclude; the pool genes carrying them are listed
    choice      the pool ranked by Spearman with the measured PSD95 punctum density
                (adult.synaptome) over the declared structures where both exist; the
                density_markers.n_markers highest make the term, the mean of their
                ranks per structure
    validation  the choice repeated on density_markers.n_halves random halves of
                those structures: how often each gene is chosen, and the Spearman of
                the chosen genes' mean with PSD95 on the other half, the agreement
                quoted. The full set's is optimistic, since the genes were chosen on
                it. Beside it, the same held-out agreement of fixed composites: the
                first proposal (FIRST_PROPOSAL), the marker panel (beyond.markers)
                and psd_pc1

Why genes when PSD95 density is measured: the measured map is one mouse and covers
about half the declared structures, while the genes are measured in nearly all of
them; PSD95 is the yardstick that chooses them, not the term. Why postsynaptic genes:
mRNA sits in somata, and a postsynaptic gene's is made by the neurons that receive
the synapses, the side the receptor sits on, whereas a presynaptic marker's marks
the neurons that send them.

The annotation is the GO Consortium's mouse GAF (MGI's) and go-basic.obo of one GO
release, downloaded once into reference/go/<release>/ under the data root and checked
against their SHA-256; the fetch log records the release and each file's own date.
Annotations with a NOT qualifier are left out; every other qualifier and every
evidence code counts. Each GO id of the rule is checked against its name in that
release, and a name that differs stops the run.

Run by run_density_markers.py.
"""

import gzip
import hashlib
import io
import re
import urllib.request
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from sepmap.adult.beyond_density import composite, first_pc
from sepmap.config import DATA, SETTINGS
from sepmap.ish import gene_sets, gene_table
from sepmap.ish.numbers import numbers_frame, numbers_path
from sepmap.ish.spatial_null import BAND
from sepmap.structures import ISH_OUT

# the rule's thresholds, the number of genes and of halves, and the genes it gave; the
# gene that encodes the stained protein and the voxels a gene value needs; the marker
# panel compared with
DENSITY_MARKERS = SETTINGS["density_markers"]
ISH = SETTINGS["ish"]
BEYOND = SETTINGS["beyond"]

OUT = ISH_OUT / "density_markers"
CANDIDATES = OUT / "candidates.csv"
EXCLUDED = OUT / "excluded.csv"
AGREEMENT = OUT / "agreement.csv"
STRUCTURES = OUT / "structures.csv"
HALVES = OUT / "halves.csv"
SELECTION = OUT / "selection.csv"
COMPARISON = OUT / "comparison.csv"
NUMBERS = numbers_path("density_markers")

# the GO release whose files are read, from the GO Consortium's archive of releases, and
# the SHA-256 of each file of it (the mouse annotation is MGI's, as GO ships it)
GO_RELEASE = "2026-08-05"
GO_ARCHIVE = f"https://release.geneontology.org/{GO_RELEASE}"
FILES = {
    "MOUSE-mod.gaf.gz": (
        "annotations/gaf",
        "6cbefce3a2efd9862f102a2732e923fde525c8fa59a3e5ce22e8178c3ebb0093",
    ),
    "go-basic.obo": (
        "ontology",
        "b08d45b268b8c24ccb2513dbbbc7d4df9f6521c099b413f79eb31e06e0fa3bcc",
    ),
}
REFERENCE = DATA / "reference" / "go" / GO_RELEASE
GAF = REFERENCE / "MOUSE-mod.gaf.gz"
OBO = REFERENCE / "go-basic.obo"
FETCH_LOG = REFERENCE / "fetch_log.txt"

CITATION = (
    "The Gene Ontology Consortium (2023). The Gene Ontology knowledgebase in 2023. "
    "Genetics 224(1), iyad031"
)

# seconds a download may take; how the downloads introduce themselves (the Gene
# Ontology's server refuses Python's default user agent)
TIMEOUT = 300
USER_AGENT = gene_table.USER_AGENT

# the cellular components that make a gene postsynaptic, with their names in GO
POOL_TERMS = {
    "GO:0014069": "postsynaptic density",
    "GO:0099572": "postsynaptic specialization",
}

# the terms of AMPA receptors and of receptor placement at the synapse: a gene annotated
# to one of them, or below it, sets how many receptors reach the surface, so it is left
# out of density
EXCLUSION_TERMS = {
    "GO:0032281": "AMPA glutamate receptor complex",
    "GO:0097113": "AMPA glutamate receptor clustering",
    "GO:2000311": "regulation of AMPA receptor activity",
    "GO:0008328": "ionotropic glutamate receptor complex",
    "GO:0035255": "ionotropic glutamate receptor binding",
    "GO:0099645": "neurotransmitter receptor localization to postsynaptic "
    "specialization membrane",
    "GO:0098696": "regulation of neurotransmitter receptor localization to "
    "postsynaptic specialization membrane",
    "GO:0099072": "regulation of postsynaptic membrane neurotransmitter receptor levels",
    "GO:0098877": "neurotransmitter receptor transport to plasma membrane",
    "GO:0099637": "neurotransmitter receptor transport",
}

# the four AMPA receptor subunits, the localisation set's panel role and the AMPA
# receptor complex family, each left out of density whatever its annotation
SUBUNITS = ("Gria1", "Gria2", "Gria3", "Gria4")
LOCALISATION_ROLE = gene_sets.PANEL_ROLE_SETS["localisation"]
AMPA_FAMILY = gene_sets.AMPA_FAMILY

# the kinds of exclusion, as the tables write them
RULE_SUBUNIT = "AMPA receptor subunit"
RULE_LOCALISATION = "localisation set"
RULE_FAMILY = "AMPA receptor complex family"
RULE_GO = "GO term"

# terms of plasticity, listed for the record and never a reason to exclude
PLASTICITY = re.compile(r"synaptic plasticity|long-term (synaptic )?potentiation", re.I)

# the composites the chosen genes are compared with: the first proposal for the term
# (two of its genes are excluded by the rule), the synaptic marker panel of [beyond], and
# the first component of the ontology panel's postsynaptic-density genes
FIRST_PROPOSAL = ("Dlg4", "Homer1", "Camk2a")
PSD_ROLE = "control_psd"
CHOSEN = "chosen"
COMPARED = {
    CHOSEN: "the genes the rule chose",
    "first_proposal": "Dlg4, Homer1, Camk2a",
    "marker_panel": "the synaptic marker panel (beyond.markers)",
    "psd_pc1": "first component of the postsynaptic-density genes",
}


# ===== Fetching =====


def sha256(data: bytes) -> str:
    """The SHA-256 of some bytes, in hex."""
    return hashlib.sha256(data).hexdigest()


def file_url(name: str) -> str:
    """Where one file of FILES lies in the GO release's archive."""
    folder, _ = FILES[name]
    return f"{GO_ARCHIVE}/{folder}/{name}"


def file_ok(path: Path, digest: str) -> bool:
    """Whether `path` exists and holds exactly the file of the release."""
    return path.exists() and sha256(path.read_bytes()) == digest


def file_version(name: str, data: bytes) -> str:
    """The version a file states in its header.

    For the GAF, its date-generated and the ontology it was made with; for the obo,
    its data-version.
    """
    if name.endswith(".gaf.gz"):
        found = []

        # the header only, read line by line from the compressed bytes
        with gzip.open(io.BytesIO(data), "rt", encoding="utf-8") as fh:
            for line in fh:
                if not line.startswith("!"):
                    break
                for key in ("date-generated:", "go-version:"):
                    if line.startswith("!" + key):
                        found.append(line[1:].strip())
        return "; ".join(found)
    for line in data.decode("utf-8").splitlines()[:5]:
        if line.startswith("data-version:"):
            return line.strip()
    return ""


def log_header() -> str:
    """The fetch log's first lines: the release, where it comes from, what to cite."""
    lines = [
        f"# Gene Ontology release {GO_RELEASE}: the mouse annotation (MGI's GAF, "
        "MOUSE-mod) and go-basic.obo",
        f"# archive     {GO_ARCHIVE}/",
        "# read by     mapping/sepmap/adult/density_markers.py, which pins the "
        "release and each file's SHA-256",
        f"# cite        {CITATION}",
        "fetched_utc\tfile\tbytes\tsha256\tversion\turl",
    ]
    return "\n".join(lines) + "\n"


def fetch(offline: bool = False) -> list[str]:
    """Download the files of FILES that are missing or differ from the release.

    A file already on disk with the release's hash is left alone; a download with
    another hash is refused and nothing is written. Each download adds a line to
    the fetch log (date, file, size, hash, the file's own version, URL). With
    `offline` a missing file stops the run instead. Returns the names downloaded.
    """
    REFERENCE.mkdir(parents=True, exist_ok=True)
    if not FETCH_LOG.exists():
        FETCH_LOG.write_text(log_header(), encoding="utf-8")
    fetched = []
    for name, (_, digest) in FILES.items():
        path = REFERENCE / name
        if file_ok(path, digest):
            continue
        if offline:
            raise FileNotFoundError(
                f"{path} is missing or is not the file of GO release {GO_RELEASE}; "
                "run without --offline to download it"
            )
        url = file_url(name)
        request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        with urllib.request.urlopen(request, timeout=TIMEOUT) as fh:
            data = fh.read()
        if sha256(data) != digest:
            raise RuntimeError(
                f"{url} has SHA-256 {sha256(data)}, not {digest}: the download was "
                "cut or the release's file changed; nothing was written"
            )
        path.write_bytes(data)
        stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%MZ")
        version = file_version(name, data)
        with open(FETCH_LOG, "a", encoding="utf-8") as fh:
            fh.write(f"{stamp}\t{name}\t{len(data)}\t{digest}\t{version}\t{url}\n")
        fetched.append(name)
    return fetched


def check_files() -> None:
    """Stop unless every file of FILES is on disk and is the file of the release."""
    for name, (_, digest) in FILES.items():
        path = REFERENCE / name
        if not file_ok(path, digest):
            raise FileNotFoundError(
                f"{path} is missing or is not the file of GO release {GO_RELEASE}: run "
                "run_density_markers.py without --offline"
            )


def versions() -> dict[str, str]:
    """{file: the version its header states}, from the files on disk."""
    return {name: file_version(name, (REFERENCE / name).read_bytes()) for name in FILES}


# ===== The annotation =====


def read_gaf(path: Path, genes: set[str]) -> tuple[pd.DataFrame, int]:
    """The annotations of `genes` in a GAF file, and how many NOT ones were left out.

    One row per annotation line of a gene of `genes` (by its symbol, column 3):
    symbol, mgi_id, qualifier, go_id, evidence, aspect and assigned_by. A line whose
    qualifier holds NOT says the gene is not there, so it is left out.
    """
    rows, n_not = [], 0
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("!"):
                continue
            field = line.rstrip("\n").split("\t")
            if field[2] not in genes:
                continue
            if "NOT" in field[3].split("|"):
                n_not += 1
                continue
            rows.append(
                dict(
                    symbol=field[2],
                    mgi_id=field[1],
                    qualifier=field[3],
                    go_id=field[4],
                    evidence=field[6],
                    aspect=field[8],
                    assigned_by=field[14],
                )
            )
    return pd.DataFrame(rows), n_not


def check_gaf(annotations: pd.DataFrame) -> None:
    """Stop when one symbol stands for two genes of the GAF (two MGI ids)."""
    ids = annotations.groupby("symbol")["mgi_id"].nunique()
    shared = sorted(ids.index[ids > 1])
    if shared:
        raise ValueError(
            f"{GAF.name}: the symbols {', '.join(shared)} stand for more than one MGI "
            "gene, so a gene of the gene table cannot be matched by its symbol"
        )


def children_of(parents: dict[str, set[str]]) -> dict[str, set[str]]:
    """GO's relations turned round: {term: the terms directly below it}."""
    children = defaultdict(set)
    for term, above in parents.items():
        for parent in above:
            children[parent].add(term)
    return dict(children)


def descendants(term: str, children: dict[str, set[str]]) -> set[str]:
    """`term` and every term below it, by the relations in `children`."""
    out, todo = set(), [term]
    while todo:
        t = todo.pop()
        if t in out:
            continue
        out.add(t)
        todo.extend(children.get(t, ()))
    return out


def check_terms(terms: dict[str, str], names: dict[str, str]) -> None:
    """Stop unless every GO id of `terms` carries its name in the ontology read.

    A rule written with an id and a name means the term of that name; when the
    release names the id otherwise, the error gives the id of the term that holds
    the name, if one does, so the rule can be corrected by hand.
    """
    for term, name in terms.items():
        if names.get(term) == name:
            continue
        holders = sorted(t for t, n in names.items() if n == name)
        raise ValueError(
            f"{OBO.name}: {term} is named {names.get(term, 'nothing')!r}, not "
            f"{name!r}; the term of that name is {', '.join(holders) or 'none'}. "
            "Correct the rule in density_markers.py and record why"
        )


def term_text(term: str, names: dict[str, str]) -> str:
    """A GO id with its name: 'GO:0014069 postsynaptic density'."""
    return f"{term} {names.get(term, '')}".strip()


def annotated_terms(annotations: pd.DataFrame) -> dict[str, dict[str, set[str]]]:
    """{gene: {GO id: the evidence codes of its annotations to it}}."""
    out = defaultdict(lambda: defaultdict(set))
    for symbol, term, evidence in zip(
        annotations["symbol"], annotations["go_id"], annotations["evidence"]
    ):
        out[symbol][term].add(evidence)
    return out


# ===== The pool =====


def exclusion_rows(
    symbol: str,
    terms: dict[str, set[str]],
    role: str,
    below: dict[str, set[str]],
    names: dict[str, str],
) -> list[dict]:
    """Every reason one gene is left out of density, one row each.

    `terms` holds the gene's GO ids with their evidence codes, `role` its role in
    the ontology panel, `below` each term of EXCLUSION_TERMS with the terms below
    it. A row: rule, term (the excluding term, or the set), via (the gene's own term
    that lies at or below it) and evidence.
    """
    rows = []
    if symbol in SUBUNITS:
        rows.append(dict(rule=RULE_SUBUNIT, term="Gria1 to Gria4", via="", evidence=""))
    if role == LOCALISATION_ROLE:
        rows.append(
            dict(
                rule=RULE_LOCALISATION,
                term=gene_sets.GENE_SETS["localisation"],
                via="",
                evidence="",
            )
        )
    if symbol in AMPA_FAMILY:
        rows.append(
            dict(
                rule=RULE_FAMILY,
                term="Schwenk et al. 2012, GO:0032281 and Gria2 to Gria4 (ish.gene_sets)",
                via="",
                evidence="",
            )
        )
    for target, covered in below.items():
        for term in sorted(set(terms) & covered):
            rows.append(
                dict(
                    rule=RULE_GO,
                    term=term_text(target, names),
                    via=term_text(term, names),
                    evidence=" ".join(sorted(terms[term])),
                )
            )
    return rows


def coverage_structures(
    profiles: dict[str, dict[str, float]], declared: list[str]
) -> list[str]:
    """The declared structures where the gene of the stained protein is measured."""
    gene = ISH["control_gene"]
    return [s for s in declared if s in profiles.get(gene, {})]


def not_in_pool_reason(row: dict, n_reference: int) -> str:
    """Why a gene is not in the pool, by the first criterion of the rule it fails."""
    if row["n_experiments_used"] < DENSITY_MARKERS["min_experiments"]:
        return (
            f"{row['n_experiments_used']} usable Allen experiments, fewer than "
            f"{DENSITY_MARKERS['min_experiments']}"
        )
    if not row["enough_coverage"]:
        return (
            f"measured in {row['n_measured']} of the {n_reference} declared structures "
            f"where Gria1 is, under {DENSITY_MARKERS['min_coverage']:.0%}"
        )
    if not row["in_gaf"]:
        return "no annotation in the mouse GAF"
    if not row["postsynaptic"]:
        return "not annotated to the postsynaptic density or specialization"
    if row["excluded"]:
        return f"excluded: {row['exclude_rules']}"
    return ""


def candidate_table(
    genes: pd.DataFrame,
    profiles: dict[str, dict[str, float]],
    declared: list[str],
    terms: dict[str, dict[str, set[str]]],
    parents: dict[str, set[str]],
    names: dict[str, str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """candidates.csv and excluded.csv: every gene of the gene table against the rule.

    `genes` is the gene table's one row per gene (ish.gene_table.per_gene),
    `profiles` every gene's merged profile, `terms` each gene's GO ids with their
    evidence (annotated_terms), `parents` and `names` the ontology. The candidates:
    a row per gene with each criterion (n_experiments_used, n_measured, coverage,
    in_gaf, postsynaptic and the terms that make it so, excluded and the rules),
    the plasticity terms it carries, eligible (every criterion but the exclusion),
    in_pool and the reason when not. The exclusions: a row per reason (exclusion_rows)
    of every gene excluded, with whether it was otherwise eligible.
    """
    children = children_of(parents)
    pool_below = set()
    for term in POOL_TERMS:
        pool_below |= descendants(term, children)
    below = {term: descendants(term, children) for term in EXCLUSION_TERMS}
    reference = coverage_structures(profiles, declared)
    role = dict(zip(genes["symbol"], genes["ontology_role"]))

    rows, excluded = [], []
    for r in genes.itertuples():
        mine = terms.get(r.symbol, {})
        measured = sum(s in profiles.get(r.symbol, {}) for s in reference)
        reasons = exclusion_rows(r.symbol, mine, role[r.symbol], below, names)
        postsynaptic = sorted(set(mine) & pool_below)
        plasticity = sorted(t for t in mine if PLASTICITY.search(names.get(t, "")))
        row = dict(
            symbol=r.symbol,
            ontology_role=r.ontology_role,
            n_experiments_used=int(r.n_experiments_used),
            n_measured=measured,
            coverage=measured / len(reference),
            enough_coverage=measured >= DENSITY_MARKERS["min_coverage"] * len(reference),
            in_gaf=bool(mine),
            postsynaptic=bool(postsynaptic),
            postsynaptic_terms="; ".join(term_text(t, names) for t in postsynaptic),
            excluded=bool(reasons),
            exclude_rules="; ".join(sorted({x["rule"] for x in reasons})),
            plasticity_terms="; ".join(term_text(t, names) for t in plasticity),
        )
        row["eligible"] = (
            row["n_experiments_used"] >= DENSITY_MARKERS["min_experiments"]
            and row["enough_coverage"]
            and row["postsynaptic"]
        )
        row["in_pool"] = row["eligible"] and not row["excluded"]
        row["reason"] = not_in_pool_reason(row, len(reference))
        rows.append(row)
        for x in reasons:
            excluded.append(dict(symbol=r.symbol, eligible=row["eligible"], **x))
    return pd.DataFrame(rows), pd.DataFrame(excluded)


# ===== The choice =====


def gene_rho(
    gene: dict[str, float], psd95: pd.Series, structures: list[str]
) -> tuple[int, float]:
    """(n, Spearman) of a gene's profile with PSD95 over the `structures` both have."""
    both = [s for s in structures if s in gene and np.isfinite(psd95.get(s, np.nan))]
    if len(both) < 3:
        return len(both), float("nan")
    rho = spearmanr([gene[s] for s in both], psd95[both].to_numpy(float)).statistic
    return len(both), float(rho)


def agreement_table(
    candidates: pd.DataFrame,
    profiles: dict[str, dict[str, float]],
    psd95: pd.Series,
    structures: list[str],
) -> pd.DataFrame:
    """agreement.csv: each eligible gene's Spearman with PSD95, the pool ranked.

    Every gene that meets the rule but for the exclusion, so the figure can mark the
    excluded ones: symbol, in_pool, excluded, exclude_rules, n_structures, rho,
    rank (within the pool, 1 the highest; empty for an excluded gene) and chosen.
    """
    eligible = candidates[candidates["eligible"]]
    rows = []
    for r in eligible.itertuples():
        n, rho = gene_rho(profiles[r.symbol], psd95, structures)
        rows.append(
            dict(
                symbol=r.symbol,
                in_pool=r.in_pool,
                excluded=r.excluded,
                exclude_rules=r.exclude_rules,
                n_structures=n,
                rho=rho,
            )
        )
    table = pd.DataFrame(rows).sort_values(["rho", "symbol"], ascending=[False, True])
    table["rank"] = np.nan
    pool = table["in_pool"].to_numpy(bool)
    table.loc[pool, "rank"] = np.arange(1, pool.sum() + 1)
    table["chosen"] = table["rank"] <= DENSITY_MARKERS["n_markers"]
    return table.reset_index(drop=True)


def chosen_genes(agreement: pd.DataFrame) -> list[str]:
    """The genes the rule chose, highest agreement first."""
    return list(agreement.loc[agreement["chosen"], "symbol"])


def check_chosen(chosen: list[str]) -> None:
    """Stop unless the rule's genes are those of settings.toml (density_markers.chosen).

    The settings hold the genes the rule gave when it was committed; a rule that now
    gives others means its inputs changed, which a reader must see before any model
    runs on them.
    """
    named = list(DENSITY_MARKERS["chosen"])
    if sorted(chosen) != sorted(named):
        raise ValueError(
            f"the rule chose {', '.join(chosen)}, but settings.toml [density_markers] "
            f"chosen names {', '.join(named)}: the gene table, the synaptome or the "
            "annotation changed. Check why before changing the settings"
        )


def structure_table(
    set_table: pd.DataFrame,
    profiles: dict[str, dict[str, float]],
    chosen: list[str],
    psd95: pd.Series,
) -> pd.DataFrame:
    """structures.csv: each declared structure, measured for the model or why not.

    The model needs Gria1 and every chosen gene measured (ish.min_voxels). Columns:
    structure, acronym, division, gria1_measured, missing_genes (the chosen genes not
    measured there), psd95_measured, in_model and reason.
    """
    gene = ISH["control_gene"]
    rows = []
    for r in set_table[set_table["in_set"]].itertuples():
        has_gene = r.structure in profiles.get(gene, {})
        missing = [g for g in chosen if r.structure not in profiles[g]]
        if not has_gene:
            reason = f"{gene} not measured"
        elif missing:
            reason = f"not measured: {' '.join(missing)}"
        else:
            reason = ""
        rows.append(
            dict(
                structure=r.structure,
                acronym=r.acronym,
                division=r.division,
                gria1_measured=has_gene,
                missing_genes=" ".join(missing),
                psd95_measured=bool(np.isfinite(psd95.get(r.structure, np.nan))),
                in_model=reason == "",
                reason=reason,
            )
        )
    return pd.DataFrame(rows)


# ===== Validation =====


def composite_on(
    genes: list[str] | tuple[str, ...],
    profiles: dict[str, dict[str, float]],
    structures: list[str],
) -> pd.Series:
    """The mean rank of `genes` over the `structures` where every one is measured."""
    both = [s for s in structures if all(s in profiles[g] for g in genes)]
    return pd.Series(composite(genes, profiles, both), index=both)


def psd_pc1_on(
    genes: pd.DataFrame, profiles: dict[str, dict[str, float]], structures: list[str]
) -> tuple[pd.Series, int]:
    """psd_pc1 over `structures`, and the number of genes it is made of.

    The first component of the ontology panel's postsynaptic-density genes (role
    control_psd) measured in every one of the structures, as analysis 4 builds it.
    """
    psd = [
        g
        for g in genes.loc[genes["ontology_role"] == PSD_ROLE, "symbol"]
        if all(s in profiles[g] for s in structures)
    ]
    pc, _ = first_pc(psd, profiles, structures)
    return pd.Series(pc, index=structures), len(psd)


def rho_with(values: pd.Series, psd95: pd.Series) -> tuple[int, float]:
    """(n, Spearman) of a composite with PSD95 over the structures both have."""
    pair = pd.DataFrame(dict(a=values, b=psd95.reindex(values.index))).dropna()
    if len(pair) < 3:
        return len(pair), float("nan")
    return len(pair), float(spearmanr(pair["a"], pair["b"]).statistic)


def fixed_composites(
    genes: pd.DataFrame,
    profiles: dict[str, dict[str, float]],
    structures: list[str],
) -> tuple[dict[str, tuple[str, ...] | pd.Series], dict[str, str]]:
    """The composites the chosen genes are compared with, and what each is made of.

    By name: the first proposal and the marker panel as gene lists, ranked within
    whatever structures they are taken on, and psd_pc1 as a map. psd_pc1 is computed
    once, over the `structures` where Gria1 is measured: it is made without PSD95, so
    a half can take it as it is.
    """
    on = [s for s in structures if s in profiles[ISH["control_gene"]]]
    pc1, n_psd = psd_pc1_on(genes, profiles, on)
    fixed = {
        "first_proposal": FIRST_PROPOSAL,
        "marker_panel": tuple(BEYOND["markers"]),
        "psd_pc1": pc1,
    }
    members = {
        "first_proposal": " ".join(FIRST_PROPOSAL),
        "marker_panel": " ".join(BEYOND["markers"]),
        "psd_pc1": f"{n_psd} genes of role {PSD_ROLE}",
    }
    return fixed, members


def composite_values(
    composite_def: tuple[str, ...] | pd.Series,
    profiles: dict[str, dict[str, float]],
    structures: list[str],
) -> pd.Series:
    """A composite on some structures: a gene list's mean rank, or a map's values."""
    if isinstance(composite_def, pd.Series):
        return composite_def.reindex([s for s in structures if s in composite_def.index])
    return composite_on(composite_def, profiles, structures)


def choose_on(
    pool: list[str],
    profiles: dict[str, dict[str, float]],
    psd95: pd.Series,
    structures: list[str],
) -> list[str]:
    """The rule's choice on some structures: the pool genes highest with PSD95 there."""
    scored = []
    for gene in pool:
        _, rho = gene_rho(profiles[gene], psd95, structures)
        scored.append((-rho, gene))
    scored.sort()
    return [gene for _, gene in scored[: DENSITY_MARKERS["n_markers"]]]


def half_split(
    pool: list[str],
    profiles: dict[str, dict[str, float]],
    psd95: pd.Series,
    structures: list[str],
    fixed: dict[str, tuple[str, ...] | pd.Series],
    seed: int = 0,
) -> pd.DataFrame:
    """halves.csv: the choice on a random half, its agreement on the other half.

    density_markers.n_halves times the `structures` are cut at random into two
    halves; the rule chooses on the first, and on the second each composite's
    Spearman with PSD95 is taken: the chosen genes' mean rank and each of `fixed`
    (fixed_composites; a gene list ranked within the half). Columns: half, genes
    (those chosen, joined), and n_ and rho_ of the chosen and of each of `fixed`.
    """
    rng = np.random.default_rng(seed)
    n_first = len(structures) // 2
    rows = []
    for half in range(DENSITY_MARKERS["n_halves"]):
        order = rng.permutation(len(structures))
        first = [structures[i] for i in sorted(order[:n_first])]
        second = [structures[i] for i in sorted(order[n_first:])]
        chosen = choose_on(pool, profiles, psd95, first)
        row = dict(half=half, genes=" ".join(chosen))
        composites = {CHOSEN: tuple(chosen), **fixed}
        for name, composite_def in composites.items():
            values = composite_values(composite_def, profiles, second)
            row[f"n_{name}"], row[f"rho_{name}"] = rho_with(values, psd95)
        rows.append(row)
    return pd.DataFrame(rows)


def selection_table(halves: pd.DataFrame, pool: list[str]) -> pd.DataFrame:
    """selection.csv: how often each pool gene is chosen on a half, most often first."""
    counts = pd.Series(0, index=pool)
    for text in halves["genes"]:
        for gene in text.split():
            counts[gene] += 1
    table = pd.DataFrame(
        dict(symbol=counts.index, n_chosen=counts.to_numpy(int))
    ).sort_values(["n_chosen", "symbol"], ascending=[False, True])
    table["share"] = table["n_chosen"] / len(halves)
    return table.reset_index(drop=True)


def comparison_table(
    chosen: list[str],
    profiles: dict[str, dict[str, float]],
    psd95: pd.Series,
    structures: list[str],
    halves: pd.DataFrame,
    fixed: dict[str, tuple[str, ...] | pd.Series],
    members: dict[str, str],
) -> pd.DataFrame:
    """comparison.csv: each composite's agreement with PSD95, held out and on the set.

    Per composite, the chosen genes first, then those of `fixed` with what each is
    made of (`members`, fixed_composites): the full set's n and Spearman (optimistic
    for the chosen genes, which were chosen on it), and over the halves the median
    held-out Spearman with its 95% range.
    """
    composites = {CHOSEN: tuple(chosen), **fixed}
    members = {CHOSEN: " ".join(chosen), **members}
    rows = []
    for name, composite_def in composites.items():
        n, rho = rho_with(composite_values(composite_def, profiles, structures), psd95)
        held = halves[f"rho_{name}"].dropna().to_numpy(float)
        lo, hi = np.percentile(held, BAND)
        rows.append(
            dict(
                composite=name,
                what=COMPARED.get(name, name),
                genes=members[name],
                n_structures=n,
                rho_full_set=rho,
                rho_held_out_median=float(np.median(held)),
                rho_held_out_lo=float(lo),
                rho_held_out_hi=float(hi),
                n_held_out=int(np.median(halves[f"n_{name}"])),
            )
        )
    return pd.DataFrame(rows)


# ===== The numbers =====


def numbers_table(
    candidates: pd.DataFrame,
    excluded: pd.DataFrame,
    agreement: pd.DataFrame,
    selection: pd.DataFrame,
    comparison: pd.DataFrame,
    counts: dict[str, int],
    releases: dict[str, str],
) -> pd.DataFrame:
    """numbers_density_markers.csv: the numbers of this step that the text quotes.

    `counts` holds the structure counts of the run (declared, reference, psd95,
    model), `releases` each GO file's own version.
    """
    eligible = candidates[candidates["eligible"]]
    pool = candidates[candidates["in_pool"]]
    out_eligible = sorted(set(excluded.loc[excluded["eligible"], "symbol"]))
    rows = [
        ("go_release", GO_RELEASE, "GO release whose GAF and go-basic.obo are read"),
        ("gaf_version", releases["MOUSE-mod.gaf.gz"], "the GAF's own header"),
        ("obo_version", releases["go-basic.obo"], "go-basic.obo's own header"),
        ("genes", len(candidates), "genes of the gene table"),
        (
            "genes_min_experiments",
            int(
                (
                    candidates["n_experiments_used"] >= DENSITY_MARKERS["min_experiments"]
                ).sum()
            ),
            "with enough usable Allen experiments",
        ),
        (
            "genes_eligible",
            len(eligible),
            "of them measured widely enough and annotated postsynaptic",
        ),
        ("genes_excluded", len(out_eligible), "of those, excluded as AMPA-linked"),
        ("excluded_genes", " ".join(out_eligible), "the genes excluded"),
        ("pool", len(pool), "genes in the pool"),
        (
            "pool_plasticity",
            int((pool["plasticity_terms"] != "").sum()),
            "pool genes carrying a plasticity term",
        ),
    ]
    rows += [(f"structures_{k}", v, f"structures: {k}") for k, v in counts.items()]
    for r in agreement[agreement["chosen"]].itertuples():
        rows.append(
            (
                f"rho_{r.symbol}",
                round(r.rho, 3),
                f"{r.symbol} with PSD95, n = {r.n_structures}",
            )
        )
    for r in comparison.itertuples():
        rows += [
            (f"{r.composite}_full_set", round(r.rho_full_set, 3), f"{r.what}, full set"),
            (
                f"{r.composite}_held_out",
                round(r.rho_held_out_median, 3),
                "held out, median",
            ),
            (
                f"{r.composite}_held_out_lo",
                round(r.rho_held_out_lo, 3),
                "2.5th percentile",
            ),
            (
                f"{r.composite}_held_out_hi",
                round(r.rho_held_out_hi, 3),
                "97.5th percentile",
            ),
        ]
    for r in selection.head(DENSITY_MARKERS["n_markers"] + 3).itertuples():
        rows.append((f"chosen_share_{r.symbol}", round(r.share, 3), "share of halves"))
    return numbers_frame(rows)
