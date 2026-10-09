"""The measured synapse density: PSD95 puncta per structure in one adult mouse (Zhu 2018).

Analysis 4 asks how much of the nano map Gria1 and synapse density leave. Its
density term is Allen mRNA of postsynaptic genes, and mRNA sits in cell bodies;
Zhu et al. 2018 counted excitatory synapses where they are. Their map is one mouse
and covers about half the declared structures, so it is the yardstick that chooses
the density genes (adult.density_markers) and a check row of analysis 4, not its
term. In a knock-in mouse with PSD95 and SAP102 tagged (Dlg4-eGFP, Dlg3-mKO2), every
punctum was detected in coronal sections of one adult male (postnatal day 80, five
18 um sections, as Hansen et al. describe it), sorted into 37 subtypes by its
intensity, size and shape, and the density of each subtype (puncta per unit area)
was measured in the regions of the Allen Reference Atlas. The table used is the one
Hansen et al. share (PLOS Biology, doi 10.1371/journal.pbio.3003637): 37 subtypes by
775 samples, a sample being one region of one hemisphere in one section, named as
'main_Thalamus_VPM left' or 'ctxall_Isocortex_layers_Layer1 left_intersect_...'.

The subtypes, in the rows of the sheet (Hansen's scpt_remap_synaptome.py indexes them
the same way): 1 to 11 hold PSD95 alone, 12 to 18 SAP102 alone, 19 to 37 both
(beyond.psd95_subtypes and beyond.sap102_subtypes). The densities taken from them,
each the mean over its subtypes:

    psd95        the 30 subtypes whose puncta hold PSD95 (1 to 11 and 19 to 37): the
                 count of PSD95 puncta, the measured excitatory-synapse density that
                 the density genes are chosen by and analysis 4's check row uses
    psd95_only   1 to 11, the "PSD95 synapses" of Hansen et al.
    sap102       the 26 subtypes holding SAP102 (12 to 37)
    all_puncta   all 37, PSD95 or SAP102

Only counts enter, never a punctum's intensity or size: how much PSD95 a synapse
holds is scaffolding, the side of surface receptor the leftover may carry, so it
cannot stand for density. As shared, each subtype's density is divided by its
largest value over the 775 samples (every row runs from exactly 0, the right locus
coeruleus with no punctum, to exactly 1), so the absolute counts and their sum over
subtypes are gone. A mean over subtypes therefore weights each subtype's map alike,
whatever its share of the puncta. Everything downstream works on ranks, so only that
weighting matters, and the variants show how far it moves the order of the
structures.

A sample is placed in a structure of the adult table through the CCF 2017 ontology:

    1  its Allen id is the one the source gives (synregions_info.pkl); a sample with
       none is looked up by the acronym in its name
    2  the id lies in the structure ('structure' term of the CCF annotation) that is
       itself or its nearest ancestor. An id above the structures (PTLp, the parent of
       VISa and VISrl; Medulla, motor related) is never spread onto the structures
       below it: its samples stay unplaced
    3  the sample's id is its unit (ACAd1, CA1slm, or the structure itself), and a
       unit's density is the mean of its samples over both hemispheres and every
       section, as the nano map pools both hemispheres
    4  a structure's density is the mean of its sampled units weighted by the voxels
       of each in the 20 um CCF annotation of the adult map, so layers count by their
       volume; a unit holding another sampled unit (PAG and its nucleus ND) weighs
       only the voxels outside it. When a unit is not drawn in the annotation (the
       piriform layers, the cerebellar granular layers), the units of that structure
       are weighted alike, and its covered share is unknown
    5  a structure none of whose parts was sampled stays missing, never filled

The source's id is followed even where its name reads otherwise (one thalamic
sample, docs/ISH_ANALYSIS.md section 4.2).

A sample whose 37 densities are all zero (the right locus coeruleus) had no punctum
detected anywhere, so it was not measured rather than empty, and is dropped. The
covered share of a structure is the share of its CCF voxels in sampled units: a
structure measured only in a thin layer is flagged by it, not dropped. Most
structures outside the cortex were sampled in one section, so the agreement of the
two hemispheres, the only check of reliability one mouse allows, says nothing about
how a structure varies from section to section.

Run by run_synaptome.py.
"""

import hashlib
import urllib.request
import xml.etree.ElementTree as ET
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from sepmap.config import DATA, SETTINGS
from sepmap.ish.numbers import numbers_frame, numbers_path
from sepmap.structures import ISH_OUT

# the subtypes of each density
BEYOND = SETTINGS["beyond"]

REFERENCE = DATA / "reference" / "synaptome"
FETCH_LOG = REFERENCE / "fetch_log.txt"
OUT = ISH_OUT / "synaptome"
SAMPLES = OUT / "samples.csv"
DENSITY = OUT / "density.csv"
COVERAGE = OUT / "coverage.csv"
AGREEMENT = OUT / "agreement.csv"
NUMBERS = numbers_path("synaptome")

# the CCF 2017 ontology, the annotation's terms, and the Allen id of each index of the
# annotation volume, as the Allen Brain Cell atlas ships them beside the annotation
ATLAS = DATA / "atlas"
TERMS = ATLAS / "parcellation_term.csv"
MEMBERSHIP = ATLAS / "parcellation_to_parcellation_term_membership.csv"
PARCELLATION = ATLAS / "parcellation.csv"

# the repository of Hansen et al. at its head of 8 October 2026; the files read are the
# same, byte for byte, at its release v1.0 (cece42b), which Zenodo archives (doi
# 10.5281/zenodo.18201390)
REPOSITORY = "netneurolab/hansen_synaptome"
COMMIT = "0399525412b6f50cfdeb5904b96da7fa8e4b507c"
FOLDER = "data/synaptome/mouse_liu2018"

# the files read, with git's hash of each at that commit (git hash-object), so a
# download is checked against the commit itself
FILES = {
    "Type_density_Ricky.xlsx": "947aec915b5f493dbd1a73be6bff78653f5d910c",
    "synregions_info.pkl": "dc1f23e891a0618936fa5d96b0102d31ccb69120",
    "type_order.npy": "aa23b02fb6e3bc9612413622444e2c81da59512e",
}
DENSITY_FILE = REFERENCE / "Type_density_Ricky.xlsx"
REGIONS_FILE = REFERENCE / "synregions_info.pkl"
TYPE_ORDER_FILE = REFERENCE / "type_order.npy"

CITATION = (
    "Zhu F, Cizeron M, Qiu Z, Benavides-Piccione R, Kopanitsa MV, Skene NG, "
    "Koniaris B, DeFelipe J, Fransen E, Komiyama NH, Grant SGN (2018) Architecture "
    "of the mouse brain synaptome. Neuron 99:781-799, doi 10.1016/j.neuron.2018.07.007",
    "Hansen JY, Luppi AI, Qiu Z, Gini S, Fulcher BD, Gozzi A, Grant SGN, Misic B. "
    "Synapse types are spatially associated with regional hemodynamics in the mouse "
    "brain. PLOS Biology, doi 10.1371/journal.pbio.3003637 (the data shared)",
)

# seconds a download may take
TIMEOUT = 120

# the subtypes of Zhu 2018, numbered 1 to 37 in the rows of the density sheet
N_SUBTYPES = 37


def subtype_rows(ranges: list[list[int]]) -> tuple[int, ...]:
    """Rows of the density sheet (from 0) of the subtypes in inclusive ranges."""
    return tuple(sorted({n - 1 for lo, hi in ranges for n in range(lo, hi + 1)}))


# the subtypes holding PSD95 and those holding SAP102, and the three kinds they make
PSD95 = subtype_rows(BEYOND["psd95_subtypes"])
SAP102 = subtype_rows(BEYOND["sap102_subtypes"])
PSD95_ONLY = tuple(sorted(set(PSD95) - set(SAP102)))
SAP102_ONLY = tuple(sorted(set(SAP102) - set(PSD95)))
BOTH = tuple(sorted(set(PSD95) & set(SAP102)))

# the densities taken per sample, each the mean over these subtypes; the first is the
# one the density genes are chosen by and analysis 4's check row uses
MEASURES = {
    "psd95": PSD95,
    "psd95_only": PSD95_ONLY,
    "sap102": SAP102,
    "all_puncta": tuple(sorted(set(PSD95) | set(SAP102))),
}
MEASURE = "psd95"

# the namespaces of an .xlsx's XML parts
XLSX = {
    "m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
    "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
}

# the parent of the ontology's root
NO_PARENT = -1

# the kinds of reason a sample is not placed in a structure
REASON_NOT_MEASURED = "not measured"
REASON_NOT_CCF = "no Allen id, and its acronym is not in the CCF ontology"
REASON_ABOVE = "above the structures, never spread onto them"
REASON_NO_STRUCTURE = "in no structure of the CCF annotation"


# ===== Fetching =====


def git_blob_sha1(data: bytes) -> str:
    """The hash git gives a file's contents (git hash-object): SHA-1 of a blob header."""
    return hashlib.sha1(f"blob {len(data)}\0".encode() + data).hexdigest()


def file_url(name: str) -> str:
    """The raw URL of one file of the source folder at the pinned commit."""
    return f"https://raw.githubusercontent.com/{REPOSITORY}/{COMMIT}/{FOLDER}/{name}"


def file_ok(path: Path, blob: str) -> bool:
    """Whether `path` exists and holds exactly the file of the pinned commit."""
    return path.exists() and git_blob_sha1(path.read_bytes()) == blob


def log_header() -> str:
    """The fetch log's first lines: the source, the commit and what to cite."""
    lines = [
        "# synapse densities of Zhu et al. 2018, as shared by Hansen et al.",
        f"# repository  https://github.com/{REPOSITORY}",
        f"# commit      {COMMIT} (files as in release v1.0, Zenodo doi "
        "10.5281/zenodo.18201390)",
        f"# folder      {FOLDER}",
    ]
    lines += [f"# cite        {c}" for c in CITATION]
    lines.append("fetched_utc\tfile\tbytes\tgit_blob_sha1\tsha256\turl")
    return "\n".join(lines) + "\n"


def fetch(offline: bool = False) -> list[str]:
    """Download the files of FILES that are missing or differ from the pinned commit.

    A file already on disk with the commit's hash is left alone; a download whose
    hash is not the commit's is refused and nothing is written. Each download adds a
    line to the fetch log (date, file, size, git and SHA-256 hashes, URL). With
    `offline` a missing file stops the run instead. Returns the names downloaded.
    """
    REFERENCE.mkdir(parents=True, exist_ok=True)
    if not FETCH_LOG.exists():
        FETCH_LOG.write_text(log_header(), encoding="utf-8")
    fetched = []
    for name, blob in FILES.items():
        path = REFERENCE / name
        if file_ok(path, blob):
            continue
        if offline:
            raise FileNotFoundError(
                f"{path} is missing or is not the file of commit {COMMIT[:7]}; "
                "run without --offline to download it"
            )
        url = file_url(name)
        with urllib.request.urlopen(url, timeout=TIMEOUT) as fh:
            data = fh.read()
        if git_blob_sha1(data) != blob:
            raise RuntimeError(
                f"{url} has git hash {git_blob_sha1(data)}, not {blob} of commit "
                f"{COMMIT[:7]}: the download was cut or the file changed; nothing "
                "was written"
            )
        path.write_bytes(data)
        stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%MZ")
        sha256 = hashlib.sha256(data).hexdigest()
        with open(FETCH_LOG, "a", encoding="utf-8") as fh:
            fh.write(f"{stamp}\t{name}\t{len(data)}\t{blob}\t{sha256}\t{url}\n")
        fetched.append(name)
    return fetched


def check_files() -> None:
    """Stop unless every file of FILES is on disk and is the file of the pinned commit."""
    for name, blob in FILES.items():
        path = REFERENCE / name
        if not file_ok(path, blob):
            raise FileNotFoundError(
                f"{path} is missing or is not the file of commit {COMMIT[:7]}: "
                "run run_synaptome.py without --offline"
            )


# ===== Reading =====


def column_index(ref: str) -> int:
    """The column of a cell reference, counted from 0: A1 -> 0, AB7 -> 27."""
    n = 0
    for ch in ref:
        if not ch.isalpha():
            break
        n = n * 26 + ord(ch.upper()) - ord("A") + 1
    return n - 1


def sheet_rows(root: ET.Element, strings: list[str]) -> list[list]:
    """The cells of one sheet's XML as rows of values: a float, a string or None."""
    m = f"{{{XLSX['m']}}}"
    cells = {}
    for row in root.iter(f"{m}row"):
        r = int(row.get("r")) - 1
        for c in row.findall("m:c", XLSX):
            kind = c.get("t", "n")
            if kind == "inlineStr":
                value = "".join(t.text or "" for t in c.iter(f"{m}t"))
            else:
                v = c.find("m:v", XLSX)
                if v is None:
                    continue
                if kind == "s":
                    value = strings[int(v.text)]
                elif kind in ("str", "e"):
                    value = v.text
                else:
                    value = float(v.text)
            cells[(r, column_index(c.get("r")))] = value
    n_rows = max((r for r, _ in cells), default=-1) + 1
    n_cols = max((c for _, c in cells), default=-1) + 1
    return [[cells.get((r, c)) for c in range(n_cols)] for r in range(n_rows)]


def read_xlsx(path: Path) -> list[list[list]]:
    """Every sheet of an .xlsx file as rows of cell values, in the workbook's order.

    venv_atlas has no spreadsheet reader, and nothing is installed for one file, so
    the file is read as what it is: a zip of XML parts, the shared strings, the
    workbook's list of sheets and one part per sheet.
    """
    with zipfile.ZipFile(path) as z:
        strings = []
        if "xl/sharedStrings.xml" in z.namelist():
            root = ET.fromstring(z.read("xl/sharedStrings.xml"))
            for si in root.findall("m:si", XLSX):
                strings.append(
                    "".join(t.text or "" for t in si.iter(f"{{{XLSX['m']}}}t"))
                )
        rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        target = {r.get("Id"): r.get("Target") for r in rels}
        book = ET.fromstring(z.read("xl/workbook.xml"))
        sheets = []
        for sheet in book.find("m:sheets", XLSX):
            part = target[sheet.get(f"{{{XLSX['r']}}}id")]

            # a target is relative to xl/, or absolute from the zip's root
            part = part.lstrip("/") if part.startswith("/") else f"xl/{part}"
            sheets.append(sheet_rows(ET.fromstring(z.read(part)), strings))
    return sheets


def check_source(
    density: np.ndarray, names: list[str], info: pd.DataFrame, type_order: np.ndarray
) -> None:
    """Stop unless the source holds what the matching and the densities rely on.

    37 subtypes by one column per sample; the sample names of the density file equal
    those of synregions_info.pkl, in order; every subtype scaled to run from exactly 0
    to exactly 1 across the samples (the scaling the docstring describes); and the
    subtype order of type_order.npy, Hansen's, keeping the three blocks PSD95 alone,
    SAP102 alone and both.
    """
    if density.shape != (N_SUBTYPES, len(names)):
        raise ValueError(
            f"{DENSITY_FILE.name}: {density.shape[0]} subtypes by {density.shape[1]} "
            f"samples, expected {N_SUBTYPES} by {len(names)} sample names"
        )
    if list(info["Region_list"]) != names:
        raise ValueError(
            f"{REGIONS_FILE.name} lists other samples than {DENSITY_FILE.name}, or in "
            "another order"
        )
    lows, highs = density.min(axis=1), density.max(axis=1)
    if not (np.all(lows == 0) and np.all(highs == 1)):
        raise ValueError(
            f"{DENSITY_FILE.name}: subtypes do not each run from 0 to 1 (lowest "
            f"{lows.max()}, highest {highs.min()}), so they are not the scaled "
            "densities this module averages"
        )
    blocks = (PSD95_ONLY, SAP102_ONLY, BOTH)
    start = 0
    for block in blocks:
        part = type_order[start : start + len(block)]
        if sorted(part.tolist()) != list(block):
            raise ValueError(
                f"{TYPE_ORDER_FILE.name} does not keep the subtypes {block[0] + 1} to "
                f"{block[-1] + 1} together: the blocks of [beyond] psd95_subtypes and "
                "sap102_subtypes are not the source's"
            )
        start += len(block)


def load_source() -> tuple[np.ndarray, pd.DataFrame]:
    """The density sheet (37 subtypes x 775 samples) and the source's table of samples.

    The first sheet holds the densities, a row per subtype numbered in column A, a
    column per sample; the second the sample names. The table of samples is
    synregions_info.pkl: name (Region_list), the acronym and Allen id Hansen et al.
    gave each sample (NaN where they found none), and their major region. A pickle
    runs code when it is read, so it is read only after its hash has matched the
    commit's.
    """
    check_files()
    first, second = read_xlsx(DENSITY_FILE)[:2]
    density = np.array([row[1:] for row in first[1:]], dtype=float)
    subtypes = [int(row[0]) for row in first[1:]]
    if subtypes != list(range(1, N_SUBTYPES + 1)):
        raise ValueError(
            f"{DENSITY_FILE.name}: subtypes {subtypes}, expected 1 to {N_SUBTYPES}"
        )
    names = [row[0] for row in second[1:]]
    info = pd.read_pickle(REGIONS_FILE)
    type_order = np.load(TYPE_ORDER_FILE, allow_pickle=False)
    check_source(density, names, info, type_order)
    return density, info.reset_index(drop=True)


# ===== The ontology =====


def load_ontology() -> pd.DataFrame:
    """The CCF 2017 ontology: Allen id, acronym, name and parent id, indexed by id.

    parcellation_term.csv also holds the Allen Brain Cell atlas's own filler terms
    (ABC-Ontology-2023), which are not Allen structures; only the AllenCCF-Ontology-2017
    terms are kept. The root's parent is NO_PARENT.
    """
    terms = pd.read_csv(TERMS, keep_default_na=False)
    terms = terms[terms["label"].str.startswith("AllenCCF-Ontology-2017")]
    parent = [
        int(p.removeprefix("MBA:")) if p else NO_PARENT
        for p in terms["parent_identifier"]
    ]
    return pd.DataFrame(
        dict(
            id=[int(i.removeprefix("MBA:")) for i in terms["identifier"]],
            acronym=list(terms["acronym"]),
            name=list(terms["name"]),
            parent=parent,
        )
    ).set_index("id")


def ancestors(i: int, parent: dict[int, int]) -> list[int]:
    """`i` and every term above it, nearest first."""
    out = []
    while i != NO_PARENT and i in parent:
        out.append(i)
        i = parent[i]
    return out


def structure_ids() -> dict[int, str]:
    """{Allen id: name} of the 'structure' terms of the CCF annotation.

    These are the structures of the adult table (volumes.per_mouse.structure_terms
    names them the same way); the Allen Brain Cell atlas's filler terms are left out.
    """
    table = pd.read_csv(MEMBERSHIP, keep_default_na=False)
    keep = (table["parcellation_term_set_name"] == "structure") & table[
        "parcellation_term_label"
    ].str.startswith("AllenCCF-Ontology-2017")
    table = table[keep].drop_duplicates("parcellation_term_label")
    ids = table["parcellation_term_label"].str.rsplit("-", n=1).str[-1].astype(int)
    return dict(zip(ids, table["parcellation_term_name"]))


def id_voxels(ann: np.ndarray, parent: dict[int, int]) -> dict[int, int]:
    """Voxels of every Allen id in a CCF annotation volume, its own and those below it.

    `ann` holds parcellation indices (the 20 um grid of the adult map,
    volumes.per_mouse.annotation_20); parcellation.csv names each index's Allen id
    ('AllenCCF-Annotation-2020-<id>'), whose voxels count for it and for every term
    above it. A term the annotation does not draw, nor anything below it, has none.
    """
    table = pd.read_csv(PARCELLATION, keep_default_na=False)
    table = table[table["label"].str.startswith("AllenCCF-Annotation")]
    id_of_index = {
        int(index): int(label.rsplit("-", 1)[-1])
        for label, index in zip(table["label"], table["parcellation_index"])
    }
    voxels = {}
    indices, counts = np.unique(ann, return_counts=True)
    for index, n in zip(indices, counts):
        # index 0 is outside the brain
        leaf = id_of_index.get(int(index))
        if leaf is None:
            continue
        for i in ancestors(leaf, parent):
            voxels[i] = voxels.get(i, 0) + int(n)
    return voxels


# ===== The samples =====


def hemisphere(name: str) -> str:
    """'left' or 'right', the last word of a sample's name (quotes and spaces ignored)."""
    side = name.strip("' ").split()[-1]
    if side not in ("left", "right"):
        raise ValueError(f"no hemisphere at the end of the sample name {name!r}")
    return side


def name_acronym(name: str) -> str:
    """The acronym a sample's name gives, as the ontology would write it.

    A region of the main list is the last field before the hemisphere, its first word
    ('main_Hippocampus_DG-mo inf left' gives DG-mo); 'CUL4,5' becomes 'CUL4, 5'. A
    cortical sample is an area intersected with a layer, in either order, and gets the
    layer's suffix: Layer2-3 is 2/3, Layer6 is 6a as the source reads it.
    """
    text = name.strip("' ")
    if "_intersect_" not in text:
        acronym = text.rsplit("_", 1)[-1].split()[0]
        return acronym.replace(",", ", ") if "," in acronym else acronym
    layer, area = "", ""
    for part in text.split("_intersect_"):
        last = part.rsplit("_", 1)[-1]
        if last.startswith("Layer"):
            layer = last.split()[0].removeprefix("Layer")
        else:
            area = last.split()[0]
    layer = {"2-3": "2/3", "6": "6a"}.get(layer, layer)
    return area + layer


def sample_table(density: np.ndarray, info: pd.DataFrame) -> pd.DataFrame:
    """One row per sample: its name, hemisphere, the source's id, and its densities.

    Columns: sample (its column in the sheet, from 0), name, hemisphere,
    source_id (NaN where the source gives none), source_acronym, name_acronym,
    measured (False when every subtype is zero), and one column per MEASURES.
    """
    table = pd.DataFrame(
        dict(
            sample=np.arange(density.shape[1]),
            name=[n.strip("' ") for n in info["Region_list"]],
            hemisphere=[hemisphere(n) for n in info["Region_list"]],
            source_id=info["ara_id"].to_numpy(dtype=float),
            source_acronym=info["acronym"].fillna("").to_numpy(dtype=object),
            name_acronym=[name_acronym(n) for n in info["Region_list"]],
            measured=~(density == 0).all(axis=0),
        )
    )
    for measure, rows in MEASURES.items():
        table[measure] = density[list(rows)].mean(axis=0)
    return table


def structures_below(i: int, children: dict[int, list[int]], structures: dict) -> int:
    """How many structures lie below the term `i` (not counting `i` itself)."""
    below, todo = 0, list(children.get(i, []))
    while todo:
        j = todo.pop()
        if j in structures:
            below += 1
        else:
            todo.extend(children.get(j, []))
    return below


def place(
    i: int,
    parent: dict[int, int],
    children: dict[int, list[int]],
    structures: dict,
    acronym: dict[int, str],
) -> tuple[int, int, str]:
    """The structure an Allen id lies in (itself or its nearest ancestor), or why none.

    Returns (structure id, -1, "") when it lies in one. Otherwise (-1, above, reason),
    `above` being the id itself or its nearest ancestor with structures below it: the
    region above the structures that the sample stands for, never spread onto them
    (PTLp1, a layer the CCF 2017 annotation no longer draws, stands for PTLp, the
    parent of VISa and VISrl). A reason is its kind, a colon, and the detail.
    """
    for a in ancestors(i, parent):
        if a in structures:
            return a, -1, ""
    for a in ancestors(i, parent):
        below = structures_below(a, children, structures)
        if below and a == i:
            return -1, a, f"{REASON_ABOVE}: {acronym[i]} holds {below}"
        if below:
            detail = f"{acronym[i]}, under {acronym[a]}, which holds {below}"
            return -1, a, f"{REASON_NO_STRUCTURE}: {detail}"
    return -1, -1, f"{REASON_NO_STRUCTURE}: {acronym[i]}"


def match_samples(
    samples: pd.DataFrame, ontology: pd.DataFrame, structures: dict[int, str]
) -> pd.DataFrame:
    """Each sample's unit and structure, by Allen id first and acronym second; or why not.

    Adds: matched_by ('id', 'acronym' or ''), unit_id and unit (the sample's own id
    and acronym), structure_id, structure, above_id (for a sample above the
    structures, the region it stands for; else -1), used, and reason (empty when
    used).
    """
    parent = ontology["parent"].to_dict()
    acronym = ontology["acronym"].to_dict()
    children = {}
    for i, p in parent.items():
        children.setdefault(p, []).append(i)
    by_acronym = {a: i for i, a in acronym.items()}

    rows = []
    for r in samples.itertuples():
        out = dict(
            matched_by="", unit_id=-1, unit="", structure_id=-1, structure="", above_id=-1
        )
        if not r.measured:
            out["reason"] = f"{REASON_NOT_MEASURED}: no punctum in any subtype"
        elif not np.isnan(r.source_id) and int(r.source_id) in parent:
            out["matched_by"], out["unit_id"] = "id", int(r.source_id)
        elif r.name_acronym in by_acronym:
            out["matched_by"], out["unit_id"] = "acronym", by_acronym[r.name_acronym]
        else:
            out["reason"] = f"{REASON_NOT_CCF}: {r.name_acronym}"

        # the unit's structure, or why it has none
        if out["unit_id"] >= 0:
            out["unit"] = acronym[out["unit_id"]]
            sid, above, why = place(out["unit_id"], parent, children, structures, acronym)
            out["structure_id"], out["above_id"], out["reason"] = sid, above, why
            if sid >= 0:
                out["structure"] = structures[sid]
        out["used"] = out["reason"] == ""
        rows.append(out)
    return pd.concat([samples, pd.DataFrame(rows, index=samples.index)], axis=1)


def placed_samples() -> tuple[pd.DataFrame, pd.DataFrame, dict[int, str], tuple]:
    """Read and check the source, and place each sample in a structure of the CCF.

    Returns the sample table (sample_table and match_samples), the CCF ontology, the
    structures of the annotation by Allen id, and the density sheet's shape
    (subtypes, samples).
    """
    density, info = load_source()
    samples = sample_table(density, info)
    ontology = load_ontology()
    structures = structure_ids()
    samples = match_samples(samples, ontology, structures)
    return samples, ontology, structures, density.shape


# ===== Per structure =====


def weighted(values: pd.Series, weights: pd.Series) -> float:
    """The mean of `values` weighted by `weights`, over the units both have."""
    w = weights.reindex(values.index)
    return float((values * w).sum() / w.sum())


def unit_weights(
    unit_ids: list[int], voxels: dict[int, int], parent: dict[int, int]
) -> tuple[pd.Series, str]:
    """The weight of each sampled unit of a structure, and 'voxels' or 'equal'.

    Units weigh by their CCF voxels. A unit that holds other sampled units (PAG,
    sampled whole and in its nucleus ND) weighs only the voxels outside them, so no
    voxel counts twice. When one of the units is not drawn in the annotation its
    voxels are unknown, and every unit of the structure weighs alike.
    """
    own = {u: float(voxels.get(u, 0)) for u in unit_ids}
    if any(v == 0 for v in own.values()):
        return pd.Series(1.0, index=unit_ids), "equal"
    weights = {}
    for u in unit_ids:
        inside = [v for v in unit_ids if v != u and u in ancestors(v, parent)]

        # only the outermost of them, so a unit nested twice is taken out once
        outermost = [
            v for v in inside if not any(x in inside for x in ancestors(v, parent)[1:])
        ]
        weights[u] = max(own[u] - sum(own[v] for v in outermost), 0.0)
    return pd.Series(weights), "voxels"


def structure_density(
    samples: pd.DataFrame, voxels: dict[int, int], parent: dict[int, int]
) -> pd.DataFrame:
    """Per structure with a sampled unit: its units, weights, coverage and densities.

    A unit's density is the mean of its samples; a structure's, the weighted mean of
    its units (unit_weights; `parent` is the ontology's, {id: parent id}).
    psd95_left and psd95_right are MEASURE from one hemisphere's samples alone.
    covered_share is the share of the structure's voxels in sampled units (1 when
    the structure itself was sampled; NaN with equal weights, when a unit's voxels
    are unknown, so a structure measured in some of its parts only is not flagged
    there). Indexed by structure name.
    """
    used = samples[samples["used"]]
    rows = []
    for (sid, structure), part in used.groupby(["structure_id", "structure"]):
        units = part.groupby("unit_id")
        weights, how = unit_weights(list(units.groups), voxels, parent)
        if how == "equal":
            covered = np.nan
        elif sid in weights.index:
            covered = 1.0
        else:
            covered = weights.sum() / voxels[sid]
        row = dict(
            structure=structure,
            n_samples=len(part),
            n_units=len(weights),
            units=" ".join(sorted(part["unit"].unique())),
            weighting=how,
            covered_share=covered,
        )
        for measure in MEASURES:
            row[measure] = weighted(units[measure].mean(), weights)
        for side in ("left", "right"):
            one = part[part["hemisphere"] == side]
            if one.empty:
                row[f"{MEASURE}_{side}"] = np.nan
            else:
                row[f"{MEASURE}_{side}"] = weighted(
                    one.groupby("unit_id")[MEASURE].mean(), weights
                )
        rows.append(row)
    return pd.DataFrame(rows).set_index("structure")


def missing_reason(
    structure_id: int, samples: pd.DataFrame, parent: dict[int, int]
) -> str:
    """Why a structure has no density: only a region above it was sampled, or nothing."""
    above = set(ancestors(structure_id, parent)[1:])
    over = samples.loc[samples["above_id"].isin(above), "unit"]
    if len(over):
        names = ", ".join(sorted(set(over)))
        return f"only a region above it was sampled ({names}), never spread onto it"
    return "no sample of it or of its parts"


def density_table(
    set_table: pd.DataFrame,
    fit: list[str],
    per_structure: pd.DataFrame,
    samples: pd.DataFrame,
    ontology: pd.DataFrame,
    structures: dict[int, str],
) -> pd.DataFrame:
    """density.csv: every structure of the adult table, its density or why it has none.

    Columns: structure, acronym, division, in_set (the declared set), in_fit (the
    structures of analysis 4's fit), measured, reason, and the columns of
    structure_density (NaN where not measured).
    """
    id_of = {name: i for i, name in structures.items()}
    parent = ontology["parent"].to_dict()
    fit = set(fit)
    rows = []
    for r in set_table.itertuples():
        measured = r.structure in per_structure.index
        if measured:
            reason = ""
        elif r.structure not in id_of:
            reason = "not a structure of the CCF ontology"
        else:
            reason = missing_reason(id_of[r.structure], samples, parent)
        rows.append(
            dict(
                structure=r.structure,
                acronym=r.acronym,
                division=r.division,
                in_set=bool(r.in_set),
                in_fit=r.structure in fit,
                measured=measured,
                reason=reason,
            )
        )
    table = pd.DataFrame(rows).set_index("structure")
    return table.join(per_structure, how="left").reset_index()


def coverage_table(density: pd.DataFrame) -> pd.DataFrame:
    """coverage.csv: per division, the structures declared, fitted and measured."""
    rows = []
    for division, part in density.groupby("division", sort=False):
        declared = part[part["in_set"]]
        fit = part[part["in_fit"]]
        if declared.empty and fit.empty:
            continue
        rows.append(
            dict(
                division=division,
                declared=len(declared),
                declared_measured=int(declared["measured"].sum()),
                fit=len(fit),
                fit_measured=int(fit["measured"].sum()),
            )
        )
    return pd.DataFrame(rows)


# ===== Agreement =====


def agreement_rows(
    density: pd.DataFrame, terms: dict[str, pd.Series], structures: list[str], label: str
) -> list[dict]:
    """Spearman of each density with each term over the measured `structures`.

    `terms` holds one value per structure (a Series indexed by name); a structure
    without a value on either side is left out of that pair.
    """
    table = density.set_index("structure")
    measured = [s for s in structures if table.loc[s, "measured"]]
    rows = []
    for measure in MEASURES:
        for term, values in terms.items():
            pair = pd.DataFrame(
                dict(a=table.loc[measured, measure], b=values.reindex(measured))
            ).dropna()
            rho = spearmanr(pair["a"], pair["b"]).statistic if len(pair) > 2 else np.nan
            rows.append(
                dict(structures=label, density=measure, term=term, n=len(pair), rho=rho)
            )
    return rows


def measure_rows(density: pd.DataFrame, structures: list[str], label: str) -> list[dict]:
    """Spearman of MEASURE with each other density over the measured `structures`.

    Each density weights the subtypes' maps differently, so these say how far the
    choice of subtypes moves the order of the structures.
    """
    table = density.set_index("structure")
    measured = table.loc[[s for s in structures if table.loc[s, "measured"]]]
    rows = []
    for other in MEASURES:
        if other == MEASURE:
            continue
        rho = spearmanr(measured[MEASURE], measured[other]).statistic
        rows.append(
            dict(structures=label, density=MEASURE, term=other, n=len(measured), rho=rho)
        )
    return rows


def hemisphere_agreement(
    density: pd.DataFrame, structures: list[str]
) -> tuple[int, float]:
    """(n, Spearman) of the left and right hemispheres' MEASURE over `structures`.

    The map is one mouse, so the agreement of its two hemispheres, sampled in the
    same sections, is the only measure of its reliability there is.
    """
    left, right = f"{MEASURE}_left", f"{MEASURE}_right"
    both = density.set_index("structure").loc[structures, [left, right]].dropna()
    if len(both) < 3:
        return len(both), float("nan")
    return len(both), float(spearmanr(both[left], both[right]).statistic)


def agreement_table(
    density: pd.DataFrame,
    fit_terms: dict[str, pd.Series],
    set_terms: dict[str, pd.Series],
    fit: list[str],
    declared: list[str],
) -> pd.DataFrame:
    """agreement.csv: each density against the maps of the fit and of the declared set.

    The densities are also set against each other.

    `fit_terms` and `set_terms` hold each map, a value per structure, on the
    structures of the fit (`fit`) and on the declared set (`declared`).
    """
    rows = agreement_rows(density, fit_terms, fit, "fit")
    rows += measure_rows(density, fit, "fit")
    rows += agreement_rows(density, set_terms, declared, "declared")
    return pd.DataFrame(rows)


def sample_numbers(samples: pd.DataFrame, density: pd.DataFrame) -> list[tuple]:
    """The numbers of the samples: how many, how many placed, and why the rest not."""
    rows = [
        (
            "samples",
            len(samples),
            "samples of the source: a region, a hemisphere, a section",
        ),
        ("samples_used", int(samples["used"].sum()), "samples placed in a structure"),
        (
            "structures_sampled",
            samples.loc[samples["used"], "structure"].nunique(),
            "CCF structures with a placed sample",
        ),
        (
            "structures_measured",
            int(density["measured"].sum()),
            "of them, structures of the adult table",
        ),
    ]
    kinds = samples.loc[~samples["used"], "reason"].str.split(":").str[0]
    for key, kind in (
        ("not_measured", REASON_NOT_MEASURED),
        ("not_ccf", REASON_NOT_CCF),
        ("above", REASON_ABOVE),
        ("no_structure", REASON_NO_STRUCTURE),
    ):
        rows.append((f"samples_{key}", int((kinds == kind).sum()), f"left out: {kind}"))
    return rows


def coverage_numbers(density: pd.DataFrame) -> list[tuple]:
    """The numbers of the coverage, of the declared set and of the fit.

    How the measured structures of the fit are measured, and the divisions missing.
    """
    fit = density[density["in_fit"]]
    declared = density[density["in_set"]]
    n_fit = int(fit["measured"].sum())
    rows = [
        ("declared", len(declared), "structures of the declared set"),
        (
            "declared_measured",
            int(declared["measured"].sum()),
            "declared structures with a measured density",
        ),
        ("fit", len(fit), "structures of analysis 4's fit"),
        ("fit_measured", n_fit, "structures of the fit with a measured density"),
        (
            "fit_share",
            round(n_fit / len(fit), 4),
            "share of the fit's structures measured",
        ),
        (
            "fit_measured_under_half",
            int((fit["covered_share"] < 0.5).sum()),
            "of them, sampled in parts making up less than half the structure",
        ),
        (
            "fit_measured_equal_weights",
            int((fit["weighting"] == "equal").sum()),
            "of them, with units weighted alike (a unit not drawn in the CCF)",
        ),
        (
            "fit_measured_two_samples",
            int((fit["n_samples"] <= 2).sum()),
            "of them, measured in at most two samples",
        ),
    ]
    missing = fit[~fit["measured"]]["division"].value_counts()
    for division, n in missing.items():
        rows.append(
            (
                f"fit_missing_{division}",
                int(n),
                f"structures of the fit in {division} missing",
            )
        )
    return rows


def numbers_table(
    samples: pd.DataFrame,
    density: pd.DataFrame,
    agreement: pd.DataFrame,
    hemispheres: tuple[int, float],
) -> pd.DataFrame:
    """numbers_synaptome.csv: the numbers of this step that the text quotes."""
    rows = sample_numbers(samples, density) + coverage_numbers(density)
    rows += [
        ("hemispheres_n", hemispheres[0], "structures sampled in both hemispheres (fit)"),
        (
            "hemispheres_rho",
            round(hemispheres[1], 4),
            "psd95, left against right hemisphere (fit)",
        ),
    ]
    psd95 = agreement[agreement["density"] == MEASURE]
    for r in psd95.itertuples():
        rows.append(
            (
                f"rho_{r.structures}_{r.term}",
                round(r.rho, 4),
                f"{r.density} against {r.term}, {r.structures} structures (n = {r.n})",
            )
        )
    return numbers_frame(rows)


def load_density() -> pd.DataFrame:
    """density.csv as run_synaptome wrote it, indexed by structure.

    The densities are NaN where a structure was not measured; the text columns
    (reason, units, weighting) are empty there instead.
    """
    if not DENSITY.exists():
        raise FileNotFoundError(f"{DENSITY} not found: run run_synaptome.py first")
    table = pd.read_csv(DENSITY)
    for column in ("reason", "units", "weighting"):
        table[column] = table[column].fillna("")
    return table.set_index("structure")
