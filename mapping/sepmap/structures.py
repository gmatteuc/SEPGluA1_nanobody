"""The declared structure set of the ISH line, the zref reference it gives, and centroids.

Every comparison of the adult map with the Allen maps runs across structures, so
which structures enter decides part of the answer. Until now each analysis took
what its inputs held: every structure of the per-mouse table, with hindbrain
structures seen in one to four adults, fibre tracts and the atlas's "...,
unassigned" leftovers. The few-mice structures do not replicate across mice
(leave-one-out rho about 0.14, against 0.93 in the forebrain), and most of the
difference between Gria1 17th and 9th in the gene ranking comes from them. So the
set is declared once, by rule, before any correlation (S1, agreed on 30 September):

    a structure enters when it is grey matter (keep_structure) and every adult
    measures it: at least region_tables.min_vox20 tissue voxels and a nano mean
    above background, in structures.min_adults of the ten adults

The same list is the reference of zref (A1). Per brain,

    v(s)     = log2(mean nano of s / mean nano of the brain's isocortex)
    zref(s)  = (v(s) - median of v over the set) / (p90 - p10 of v over the set)

where v is young_vs_adult.region_plot's cref (the isocortex mean voxel-weighted over
the brain's isocortex structures). The young-against-adult tables take the median
and the spread over the structures all 17 brains share, so a new brain there moves
every adult's zref; here they depend only on the brain and the declared list. zref
is computed for every structure, in the set or not.

Centroids, for the spatial null and the gradient control of adult.beyond_controls,
are taken in one hemisphere. Averaged over both, every bilateral structure sits on
the midline and the medio-lateral distance between structures vanishes. The CCF
annotation is symmetric, so either half gives the same distances; the half with ML
index below the midline is used.

Run by run_structure_set.py.
"""

from collections.abc import Iterable
from pathlib import Path

import numpy as np
import pandas as pd

from sepmap.config import DATA, SETTINGS
from sepmap.volumes.per_mouse import annotation_20, structure_terms

# the adults a structure must be measured in, the grey-matter divisions and P9's
# divisions; the smallest structure a brain's table keeps
STRUCTURES = SETTINGS["structures"]
REGION_TABLES = SETTINGS["region_tables"]

TABLES = DATA / "adult_v2" / "ish_analysis" / "tables"
STRUCTURE_SET = TABLES / "structure_set.csv"
CENTROIDS = TABLES / "centroids.csv"

GREY = set(STRUCTURES["grey"])

# a 20 um voxel of the CCF grid, in mm
VOXEL_MM = 0.02

# the reasons a structure is left out, in the order the rule applies them
REASON_FEW_ADULTS = "measured in fewer than {} adults"
REASON_NOT_GREY = "not grey matter"
REASON_CATCH_ALL = "catch-all label, not a structure"


def keep_structure(name: str, division: str) -> tuple[bool, str]:
    """Whether a structure is grey matter and a real structure: (keep, reason why not).

    Two rules, both applied before any fitting so neither can be tuned to the
    answer: grey matter only, by division; and no "..., unassigned" entries, which
    are voxels the atlas could not place rather than structures.
    """
    if "unassigned" in name.lower():
        return False, REASON_CATCH_ALL
    if division not in GREY:
        return False, f"division {division} is {REASON_NOT_GREY}"
    return True, ""


def structure_set(per_mouse: pd.DataFrame, adults: list[str]) -> pd.DataFrame:
    """One row per structure of the adult table: in the declared set or not, and why.

    `per_mouse` holds a row per adult and structure (mouse, structure, acronym,
    division, nano_mean), for structures with at least region_tables.min_vox20
    tissue voxels in that brain. A structure counts as measured in an adult when it
    has a row there and its nano mean is above background (its log2 exists).
    Columns: structure, acronym, division, n_adults, grey, in_set, reason.
    """
    min_adults = STRUCTURES["min_adults"]
    measured = per_mouse[per_mouse["mouse"].isin(adults) & (per_mouse["nano_mean"] > 0)]
    first = per_mouse.drop_duplicates("structure").set_index("structure")
    n_adults = measured.groupby("structure")["mouse"].nunique()
    rows = []
    for structure in sorted(first.index):
        division = first.loc[structure, "division"]
        n = int(n_adults.get(structure, 0))
        grey, why_not = keep_structure(structure, division)

        # the rule in the order of the funnel: every adult first, then grey matter
        if n < min_adults:
            reason = REASON_FEW_ADULTS.format(min_adults)
        elif not grey:
            reason = why_not
        else:
            reason = ""
        rows.append(
            dict(
                structure=structure,
                acronym=first.loc[structure, "acronym"],
                division=division,
                n_adults=n,
                grey=grey,
                in_set=reason == "",
                reason=reason,
            )
        )
    return pd.DataFrame(rows)


def load_structure_set(path: Path | None = None) -> pd.DataFrame:
    """The table run_structure_set wrote, with in_set and grey as booleans."""
    if path is None:
        path = STRUCTURE_SET
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found: run run_structure_set.py first, it declares the set"
        )
    table = pd.read_csv(path, keep_default_na=False)
    table["in_set"] = table["in_set"].astype(str) == "True"
    table["grey"] = table["grey"].astype(str) == "True"
    return table


def declared_structures(path: Path | None = None) -> list[str]:
    """The structures of the declared set, sorted by name."""
    table = load_structure_set(path)
    return sorted(table.loc[table["in_set"], "structure"])


def zref(v: pd.Series, reference: Iterable[str]) -> tuple[pd.Series, float, float]:
    """zref of one brain from its v (log2 over its isocortex mean), indexed by structure.

    The median and the p90 - p10 spread are taken over the structures of
    `reference` that have a value in `v`; the spread is floored at 1e-6, as in
    young_vs_adult.region_plot.range_match, so a flat brain cannot divide by zero.
    Returns (zref of every structure of `v`, median, spread).
    """
    ref = v.reindex(sorted(set(reference))).dropna()
    if ref.empty:
        raise ValueError("no structure of the reference set has a value in this brain")
    p10, median, p90 = np.percentile(ref.to_numpy(), [10, 50, 90])
    spread = max(p90 - p10, 1e-6)
    return (v - median) / spread, float(median), float(spread)


def name_volume(ann: np.ndarray, names: dict[int, str]) -> tuple[np.ndarray, list[str]]:
    """Each voxel's structure, as an index into the returned list of structure names.

    The layers of an area are separate parcellation indices with one structure
    name, so they become one label. Code 0 is outside the brain; an index with no
    name in the ontology keeps one of its own, "id<index>", as region_plot does.
    """
    label_names = [""]
    code_of_name = {}
    lut = np.zeros(int(ann.max()) + 1, dtype=np.int32)
    for idx in np.unique(ann):
        if idx == 0:
            continue
        name = names.get(int(idx), f"id{idx}")
        if name not in code_of_name:
            code_of_name[name] = len(label_names)
            label_names.append(name)
        lut[idx] = code_of_name[name]
    return lut[ann], label_names


def border_voxels(labels: np.ndarray) -> np.ndarray:
    """Voxels with a six-neighbour of another label (outside the brain included)."""
    border = np.zeros(labels.shape, dtype=bool)
    for axis in range(labels.ndim):
        n = labels.shape[axis]
        lo = [slice(None)] * labels.ndim
        hi = [slice(None)] * labels.ndim
        lo[axis] = slice(0, n - 1)
        hi[axis] = slice(1, n)
        differs = labels[tuple(lo)] != labels[tuple(hi)]
        border[tuple(lo)] |= differs
        border[tuple(hi)] |= differs
    return border & (labels > 0)


def centroids(
    labels: np.ndarray,
    label_names: list[str],
    structures: Iterable[str],
    one_hemisphere: bool = True,
) -> pd.DataFrame:
    """Mean (AP, DV, ML) position of each structure in mm, and its voxels counted.

    `labels` is a (AP, DV, ML) volume of codes into `label_names` (name_volume)
    on the 20 um grid. With `one_hemisphere` only voxels with an ML index below
    half the grid count, so a bilateral structure sits in its own hemisphere and
    not on the midline. A structure with no voxel gets NaN.
    """
    if one_hemisphere:
        half = labels.shape[2] // 2
        labels = labels[:, :, :half]
    # the voxel count and the summed coordinates of every label, in one pass
    ap, dv, ml = np.nonzero(labels)
    codes = labels[ap, dv, ml]
    n_labels = len(label_names)
    count = np.bincount(codes, minlength=n_labels)
    sums = [np.bincount(codes, weights=c, minlength=n_labels) for c in (ap, dv, ml)]
    code_of = {name: i for i, name in enumerate(label_names) if i > 0}

    rows = []
    for structure in sorted(set(structures)):
        code = code_of.get(structure)
        if code is None or count[code] == 0:
            position = np.full(3, np.nan)
            n = 0
        else:
            n = int(count[code])
            position = np.array([s[code] for s in sums]) / n * VOXEL_MM
        rows.append(
            dict(
                structure=structure,
                ap_mm=position[0],
                dv_mm=position[1],
                ml_mm=position[2],
                n_voxels=n,
            )
        )
    return pd.DataFrame(rows)


def ccf_labels() -> tuple[np.ndarray, list[str]]:
    """The adult CCF annotation on the 20 um grid as structure codes, and their names."""
    names, _, _ = structure_terms()
    return name_volume(annotation_20("ccf"), names)


def ccf_centroids(
    structures: Iterable[str], one_hemisphere: bool = True
) -> dict[str, np.ndarray]:
    """{structure: (AP, DV, ML) in mm} on the adult CCF, one hemisphere by default."""
    labels, label_names = ccf_labels()
    table = centroids(labels, label_names, structures, one_hemisphere)
    return {
        r["structure"]: np.array([r["ap_mm"], r["dv_mm"], r["ml_mm"]])
        for _, r in table.iterrows()
    }


def load_centroids(path: Path | None = None) -> pd.DataFrame:
    """The centroids run_structure_set wrote, indexed by structure."""
    if path is None:
        path = CENTROIDS
    if not path.exists():
        raise FileNotFoundError(f"{path} not found: run run_structure_set.py first")
    return pd.read_csv(path).set_index("structure")
