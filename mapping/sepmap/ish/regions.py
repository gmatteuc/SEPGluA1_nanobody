"""Allen ISH expression grids to one mean per gene and structure.

The gene side is put on the footing of the nano tables: the same structure names,
the same ontology, one row per gene and structure, so the comparison downstream is
a join. An experiment is a file pair in atlas_ish/, `<id>_energy.mhd` (a text
header) and `<id>_energy.raw` (float32): expression energy (intensity x density)
on a 67 x 41 x 58 grid at 200 um, where -1 marks a voxel with no data (about half
of them, since the box is larger than the brain). Those are dropped, never read as
zero. A few experiments come in a box of their own (68 x 40 x 50, 73 x 41 x 53)
with the same `Offset = 0 0 0`, so nothing places them against the atlas; cropping
them would put every voxel in the wrong structure, so they are dropped and listed,
with the reason, in the drops table beside the region table.

MetaImage stores x fastest, so numpy reads (z, y, x) and `transpose(2, 1, 0)`
gives (AP, DV, ML). The orientation was checked once, outside this code: of the 48
axis permutations and flips, this one matches the grid's valid-data mask to the
atlas brain mask best (Dice 0.81, the next distinct one 0.69), and the old MATLAB
route declares the same layout. That test cannot resolve a left-right flip, which
does not matter here: a structure mean pools both hemispheres.

Voxels are assigned with the CCF annotation sampled onto the same 200 um grid, and
structures are keyed by name, summing their layer-level indices, as
young_vs_adult.region_plot does for the nano side. Two means per structure:

    full      every valid voxel inside the structure
    eroded    the same after peeling one 200 um voxel off the border

ISH is coarse and registration imperfect, so a small structure's mean is partly
its neighbours'; the eroded mean guards against that, and both are written so the
comparison can show whether it matters. A coverage column tells a structure
measured from three voxels from one measured from three thousand.

Run by run_ish_regions.py.
"""

import csv
import time
from collections import defaultdict
from pathlib import Path

import nibabel as nib
import numpy as np
from scipy.ndimage import binary_erosion

from sepmap.config import DATA, SETTINGS
from sepmap.volumes.per_mouse import structure_terms

# the valid voxels a structure needs to get a value
ISH_REGIONS = SETTINGS["ish_regions"]

ISH_DIR = DATA / "atlas_ish"
OUT = DATA / "adult_v2" / "ish"

# the panel passes of settings.toml ([ish_panels]): each names a panel and the
# table it writes, and is chosen by name (run_ish_regions.py --panel), so both
# panels go through this one aggregation. A relative panel path is taken inside the
# data root, so the same setting works on a copy of the data; an absolute one is
# used as it is
ISH_PANELS = SETTINGS["ish_panels"]

# the default pass: the original 100-gene panel
DEFAULT_PANEL = "targets"

# voxel size of the Allen grids, in um
GRID_UM = 200

# the Allen flag for "no data here"
MISSING = -1.0

# the shared reference box, in the header's (x, y, z) order
GRID_DIMS = (67, 41, 58)


def panel_files(name: str) -> tuple[Path, str]:
    """(panel CSV path, output table name) of one panel pass in settings.toml."""
    if name not in ISH_PANELS:
        raise ValueError(
            f"no ISH panel pass {name!r} in settings.toml; the passes "
            f"are {', '.join(ISH_PANELS)}"
        )
    panel = Path(ISH_PANELS[name]["panel"])
    if not panel.is_absolute():
        panel = DATA / panel
    return panel, ISH_PANELS[name]["table"]


class NotReferenceGrid(Exception):
    """The experiment was gridded in a box of its own, so it cannot be placed."""


def read_energy(experiment_id: int | str) -> np.ndarray:
    """One Allen grid as (AP, DV, ML) at 200 um, with missing voxels as NaN.

    The header is read rather than trusted: if a future download has a different
    size or spacing, this raises instead of quietly reshaping into nonsense.
    """
    mhd = ISH_DIR / f"{experiment_id}_energy.mhd"
    raw = ISH_DIR / f"{experiment_id}_energy.raw"
    hdr = {}
    with open(mhd) as fh:
        for line in fh:
            if "=" in line:
                k, v = line.split("=", 1)
                hdr[k.strip()] = v.strip()

    # the header's (x, y, z) is (AP, DV, ML)
    dims = [int(x) for x in hdr["DimSize"].split()]
    spacing = {float(x) for x in hdr["ElementSpacing"].split()}
    if spacing != {float(GRID_UM)}:
        raise ValueError(f"{experiment_id}: spacing {spacing} um, expected {GRID_UM}")
    if tuple(dims) != GRID_DIMS:
        raise NotReferenceGrid(f"grid is {tuple(dims)}, not {GRID_DIMS}")
    vol = np.fromfile(raw, dtype=np.float32)
    if vol.size != np.prod(dims):
        raise ValueError(
            f"{experiment_id}: {vol.size} values, header says {np.prod(dims)}"
        )

    # x is stored fastest, so numpy reads (z, y, x); back to (AP, DV, ML)
    vol = vol.reshape(dims[::-1]).transpose(2, 1, 0)
    return np.where(vol == MISSING, np.nan, vol)


def annotation_200() -> np.ndarray:
    """CCF labels on the ISH grid: the 10 um annotation sampled every 20th voxel.

    The result is one voxel smaller than the Allen grid in each axis (66 x 40 x 57
    against 67 x 41 x 58) because their box is slightly larger; the offset that
    aligns them is zero, which is what the orientation check measured.
    """
    ann = np.asarray(nib.load(DATA / "atlas" / "annotation_10.nii.gz").dataobj)
    return ann[::20, ::20, ::20]


def region_means(
    vol: np.ndarray, ann: np.ndarray, names: dict[int, str], eroded_ann: np.ndarray
) -> dict[str, tuple[float, float, int, int]]:
    """Full and eroded mean, valid voxels and all voxels of each structure.

    Returns {structure name: (full mean, eroded mean, n valid voxels, n voxels in
    the structure)}, summed over the layer-level indices that share a structure
    name, as on the nano side, so the two tables join on the same key. A structure
    with fewer than ish_regions.min_voxels valid voxels is left out; its eroded
    mean is NaN
    when erosion leaves fewer than that.
    """
    out = {}

    # per structure name: sum, n, eroded sum, eroded n, n voxels in total
    acc = defaultdict(lambda: [0.0, 0, 0.0, 0, 0])
    valid = np.isfinite(vol)
    for idx in np.unique(ann):
        if idx == 0 or idx not in names:
            continue
        m = ann == idx
        mv = m & valid
        a = acc[names[idx]]
        a[0] += float(vol[mv].sum())
        a[1] += int(mv.sum())
        a[4] += int(m.sum())
        me = (eroded_ann == idx) & valid
        a[2] += float(vol[me].sum())
        a[3] += int(me.sum())
    min_voxels = ISH_REGIONS["min_voxels"]
    for name, (s, n, se, ne, ntot) in acc.items():
        if n >= min_voxels:
            out[name] = (s / n, (se / ne) if ne >= min_voxels else np.nan, n, ntot)
    return out


def main(only: list[str] | None = None, panel_name: str = DEFAULT_PANEL) -> None:
    """Write the region table of one panel pass, and the table of dropped genes.

    With `only`, a list of gene symbols, only those genes of the panel.
    """
    # the panel's experiments, the structure names and the annotation on the grid
    panel_path, table_name = panel_files(panel_name)
    OUT.mkdir(parents=True, exist_ok=True)
    panel = [r for r in csv.DictReader(open(panel_path, newline="", encoding="utf-8"))]
    if only:
        want = {g.lower() for g in only}
        panel = [r for r in panel if r["symbol"].lower() in want]
    names, _, _ = structure_terms()
    ann_full = annotation_200()

    # erode each structure by one voxel, once, and reuse it for every gene
    print("eroding structure masks once (200 um, one voxel)...", flush=True)
    eroded = np.zeros_like(ann_full)
    for idx in np.unique(ann_full):
        if idx == 0:
            continue
        m = ann_full == idx
        e = binary_erosion(m)

        # a structure too small to erode keeps its full mask
        eroded[e if e.any() else m] = idx

    # one mean per structure for every experiment of the panel
    rows, dropped = [], []
    for i, gene in enumerate(panel, 1):
        t0 = time.time()
        sym, eid = gene["symbol"], gene["experiment_id"]
        try:
            vol = read_energy(eid)
        except FileNotFoundError:
            dropped.append(
                dict(symbol=sym, experiment_id=eid, reason="grid not downloaded")
            )
            print(
                f"{i:3d}/{len(panel)} {sym:10s} DROPPED  grid not downloaded", flush=True
            )
            continue
        except NotReferenceGrid as why:
            dropped.append(dict(symbol=sym, experiment_id=eid, reason=str(why)))
            print(f"{i:3d}/{len(panel)} {sym:10s} DROPPED  {why}", flush=True)
            continue
        vol = vol[: ann_full.shape[0], : ann_full.shape[1], : ann_full.shape[2]]
        rm = region_means(vol, ann_full, names, eroded)
        for name, (full, ero, n, ntot) in rm.items():
            rows.append(
                dict(
                    symbol=sym,
                    experiment_id=eid,
                    category=gene["category"],
                    plane=gene.get("plane", ""),
                    structure=name,
                    ish_mean=f"{full:.6g}",
                    ish_mean_eroded=("" if np.isnan(ero) else f"{ero:.6g}"),
                    n_voxels=n,
                    n_voxels_structure=ntot,
                    coverage=f"{n / max(ntot, 1):.3f}",
                )
            )
        print(
            f"{i:3d}/{len(panel)} {sym:10s} {len(rm):3d} structures   "
            f"{time.time() - t0:.1f} s",
            flush=True,
        )

    # the region table, and the experiments dropped with their reason
    path = OUT / table_name
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    n_genes = len({r["symbol"] for r in rows})
    print(
        f"\n{len(rows):,} rows, {len(panel) - len(dropped)} experiments, "
        f"{n_genes} genes -> {path}",
        flush=True,
    )

    path = OUT / (table_name.replace(".csv", "") + "_drops.csv")
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["symbol", "experiment_id", "reason"])
        w.writeheader()
        w.writerows(dropped)
    print(f"{len(dropped)} genes dropped -> {path}", flush=True)
