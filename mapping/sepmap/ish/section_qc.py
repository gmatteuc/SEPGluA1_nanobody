"""Section quality of every Allen ISH experiment: failed sections flagged, never filled.

An Allen expression grid is built from one series of sections of one P56 mouse,
and a section can fail: torn, out of focus, or stained far more weakly than its
neighbours. Its voxels then read as no expression where there is some, and a
structure mean that pools them is too low. P9 found such sections by comparing
each section's median with the brain-wide median section (below 0.3 of it) and
filled them by interpolating from the good sections on either side. Both halves
were wrong in places: a brain-wide threshold flags true absence (Slc17a6 is not
expressed in the anterior forebrain, and its sections 15 to 18 were "repaired"),
and interpolation invents values. So here (A2):

    per section    along the experiment's own section axis (AP for a coronal
                   series, ML for a sagittal one), inside the CCF brain sampled
                   on the 200 um grid: the in-brain voxels, the voxels with data,
                   and the median energy of those
    judged         a section with at least ish_qc.min_plane_voxels in-brain voxels
                   and some data; the first and last sections hold a sliver of
                   brain whose median is noise
    local ref      the median of the medians of the judged sections within
                   ish_qc.neighbours on each side (600 um each way): wider than one
                   failed section, narrower than the brain's large gradients
    flagged        a judged section whose median is below ish_qc.local_fraction of
                   its local reference, five-fold dimmer than its neighbours

A flagged section is set missing (NaN), as the grid's own -1 voxels are, and every
other voxel is left as it was: nothing is interpolated. The structure means then
come from the sections that remain.

Whether a flagged section is a failure or a true absence of expression is a human
call. mapping/ish_section_exceptions.csv lists the sections reviewed as true
absence (symbol, experiment_id, sections, reason, status): with status proposed or
accepted they are kept, with rejected they are set missing like any other flag.
The build proposes candidates, which stay "proposed" until reviewed on the QC
sheets (figures/qc/).

Run by run_ish_section_qc.py.
"""

import csv
from pathlib import Path

import numpy as np
import pandas as pd

from sepmap.config import SETTINGS, code_root
from sepmap.ish.regions import NotReferenceGrid, read_energy
from sepmap.structures import TABLES

# the neighbours of the local reference, the share below which a section is flagged,
# the in-brain voxels a section needs to be judged
ISH_QC = SETTINGS["ish_qc"]

SECTION_QC = TABLES / "section_qc.csv"
EXPERIMENT_QC = TABLES / "experiment_qc.csv"
EXCEPTIONS = Path(code_root()) / "mapping" / "ish_section_exceptions.csv"

# the axis of the (AP, DV, ML) grid along which each plane of section was cut
SECTION_AXIS = {"coronal": 0, "sagittal": 2}
AXIS_NAME = {0: "AP", 2: "ML"}

# exception statuses that keep a flagged section as a true absence
KEEP_STATUSES = ("proposed", "accepted")

# what a section's status says, as the tables and the sheets write it
OK = "ok"
FLAGGED = "flagged"
ABSENCE = "absence kept"
NO_DATA = "no data"
TOO_LITTLE_BRAIN = "too little brain"


def section_profile(vol: np.ndarray, brain: np.ndarray, axis: int) -> pd.DataFrame:
    """One row per section along `axis`: in-brain voxels, voxels with data, median.

    `vol` is the grid as (AP, DV, ML) with missing voxels NaN, cut to the shape of
    `brain`, the CCF brain mask on the same grid. The median is over the in-brain
    voxels with data, NaN when there are none.
    """
    rows = []
    for k in range(vol.shape[axis]):
        plane = np.take(vol, k, axis=axis)
        inside = np.take(brain, k, axis=axis)
        values = plane[inside & np.isfinite(plane)]
        n_brain = int(inside.sum())
        rows.append(
            dict(
                section=k,
                n_brain=n_brain,
                n_valid=int(values.size),
                valid_fraction=values.size / n_brain if n_brain else np.nan,
                median_energy=float(np.median(values)) if values.size else np.nan,
            )
        )
    return pd.DataFrame(rows)


def local_reference(medians: np.ndarray, judged: np.ndarray) -> np.ndarray:
    """Each section's local reference: the median of its judged neighbours' medians.

    The neighbours are the sections within ish_qc.neighbours on each side, the
    section itself left out; NaN when none of them is judged.
    """
    k_max = ISH_QC["neighbours"]
    n = len(medians)
    out = np.full(n, np.nan)
    for k in range(n):
        lo, hi = max(0, k - k_max), min(n, k + k_max + 1)
        near = [medians[j] for j in range(lo, hi) if j != k and judged[j]]
        if near:
            out[k] = float(np.median(near))
    return out


def flag_sections(profile: pd.DataFrame, kept_absent: set[int]) -> pd.DataFrame:
    """The profile with its local reference and each section's status and reason.

    `kept_absent` holds the sections an exception keeps as true absence. A status
    is ok, flagged, absence kept, no data (no voxel with data in the brain) or too
    little brain (fewer than ish_qc.min_plane_voxels in-brain voxels).
    """
    out = profile.copy()
    big = out["n_brain"].to_numpy() >= ISH_QC["min_plane_voxels"]
    has_data = out["n_valid"].to_numpy() > 0
    judged = big & has_data
    medians = out["median_energy"].to_numpy()
    ref = local_reference(medians, judged)
    out["local_reference"] = ref

    # a judged section five-fold dimmer than its neighbours, unless reviewed as absence
    dim = judged & np.isfinite(ref) & (medians < ISH_QC["local_fraction"] * ref)
    status, reason = [], []
    for k in range(len(out)):
        if not big[k]:
            status.append(TOO_LITTLE_BRAIN)
            reason.append(f"{out['n_brain'].iloc[k]} in-brain voxels")
        elif not has_data[k]:
            status.append(NO_DATA)
            reason.append("no voxel with data in the brain")
        elif dim[k] and k in kept_absent:
            status.append(ABSENCE)
            reason.append("dim, but listed as a true absence")
        elif dim[k]:
            status.append(FLAGGED)
            reason.append(
                f"median {medians[k]:.3g} below {ISH_QC['local_fraction']} x {ref[k]:.3g}"
            )
        else:
            status.append(OK)
            reason.append("")
    out["status"] = status
    out["reason"] = reason
    return out


def apply_flags(vol: np.ndarray, sections: list[int], axis: int) -> np.ndarray:
    """A copy of `vol` with the given sections missing (NaN); nothing else changes."""
    out = vol.copy()
    for k in sections:
        index = [slice(None)] * vol.ndim
        index[axis] = k
        out[tuple(index)] = np.nan
    return out


def load_exceptions(path: Path | None = None) -> dict[tuple[str, str], dict]:
    """The reviewed true absences, {(symbol, experiment id): {sections, status, reason}}.

    Only exceptions whose status keeps a section (proposed, accepted) give
    sections; a rejected one is listed with no sections, so its flags stand.
    """
    if path is None:
        path = EXCEPTIONS
    out = {}
    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r["status"] not in KEEP_STATUSES + ("rejected",):
                raise ValueError(
                    f"{path}: status {r['status']!r} for {r['symbol']} "
                    f"{r['experiment_id']}; it must be proposed, accepted or rejected"
                )
            sections = {int(s) for s in r["sections"].split()}
            if r["status"] == "rejected":
                sections = set()
            out[(r["symbol"], str(r["experiment_id"]))] = dict(
                sections=sections, status=r["status"], reason=r["reason"]
            )
    return out


def read_grid(experiment_id: str, shape: tuple[int, int, int]) -> np.ndarray:
    """One experiment's grid as (AP, DV, ML), missing voxels NaN, cut to `shape`.

    The Allen box is one voxel larger than the CCF sampled at 200 um in each axis
    (ish.regions.annotation_200), and the offset that aligns them is zero.
    """
    vol = read_energy(experiment_id)
    return vol[: shape[0], : shape[1], : shape[2]]


def experiment_qc(
    row: dict, brain: np.ndarray, exceptions: dict
) -> tuple[pd.DataFrame | None, dict]:
    """The section table of one experiment and its one-line summary.

    `row` holds symbol, experiment_id and plane. When the grid is missing or in a
    box of its own the section table is None and the summary says why.
    """
    symbol, eid, plane = row["symbol"], str(row["experiment_id"]), row["plane"]
    summary = dict(
        symbol=symbol,
        experiment_id=eid,
        plane=plane,
        grid="ok",
        axis="",
        n_judged=0,
        n_flagged=0,
        flagged_sections="",
        n_absence_kept=0,
        absence_sections="",
        exception_status="",
    )
    try:
        vol = read_grid(eid, brain.shape)
    except FileNotFoundError:
        summary.update(grid="missing")
        return None, summary
    except NotReferenceGrid as why:
        summary.update(grid=str(why))
        return None, summary
    axis = SECTION_AXIS[plane]
    exception = exceptions.get((symbol, eid), {})
    table = flag_sections(
        section_profile(vol, brain, axis), exception.get("sections", set())
    )
    table.insert(0, "axis", AXIS_NAME[axis])
    table.insert(0, "plane", plane)
    table.insert(0, "experiment_id", eid)
    table.insert(0, "symbol", symbol)
    summary.update(
        axis=AXIS_NAME[axis],
        n_judged=int(table["status"].isin([OK, FLAGGED, ABSENCE]).sum()),
        n_flagged=int((table["status"] == FLAGGED).sum()),
        flagged_sections=" ".join(
            str(k) for k in table.loc[table["status"] == FLAGGED, "section"]
        ),
        n_absence_kept=int((table["status"] == ABSENCE).sum()),
        absence_sections=" ".join(
            str(k) for k in table.loc[table["status"] == ABSENCE, "section"]
        ),
        exception_status=exception.get("status", ""),
    )
    return table, summary


def qc_tables(
    experiments: pd.DataFrame, brain: np.ndarray, exceptions: dict
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """The section table of every experiment, and one summary row per experiment."""
    tables, summaries = [], []
    for i, row in enumerate(experiments.to_dict("records"), 1):
        table, summary = experiment_qc(row, brain, exceptions)
        summaries.append(summary)
        if table is not None:
            tables.append(table)
        if i % 100 == 0:
            print(f"  {i}/{len(experiments)} experiments", flush=True)
    return pd.concat(tables, ignore_index=True), pd.DataFrame(summaries)


def flagged_sections(summary: pd.DataFrame) -> dict[str, list[int]]:
    """{experiment id: sections to set missing}, from the experiment summary."""
    out = {}
    for _, r in summary.iterrows():
        text = r.get("flagged_sections", "")
        if isinstance(text, str) and text:
            out[str(r["experiment_id"])] = [int(k) for k in text.split()]
    return out


def load_experiment_qc(path: Path | None = None) -> pd.DataFrame:
    """The experiment summary run_ish_section_qc wrote."""
    if path is None:
        path = EXPERIMENT_QC
    if not path.exists():
        raise FileNotFoundError(f"{path} not found: run run_ish_section_qc.py first")
    return pd.read_csv(path, dtype={"experiment_id": str}, keep_default_na=False)


def load_section_qc(path: Path | None = None) -> pd.DataFrame:
    """The section table run_ish_section_qc wrote."""
    if path is None:
        path = SECTION_QC
    if not path.exists():
        raise FileNotFoundError(f"{path} not found: run run_ish_section_qc.py first")
    return pd.read_csv(
        path, dtype={"experiment_id": str}, keep_default_na=False, na_values=[""]
    )
