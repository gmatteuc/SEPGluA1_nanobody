"""
Allen ISH expression grids -> one mean per gene per structure.

The point of this step is to put gene expression on exactly the same footing as
the nano measurements: the same structure names, the same ontology, one row per
(gene, structure), so the comparison downstream is a join and nothing more.

What the Allen grids are. One file pair per experiment in data\\atlas_ish:
`<id>_energy.mhd` (a text header) and `<id>_energy.raw` (float32). The grid is
67 x 41 x 58 at 200 um, values are expression ENERGY (intensity x density), and
**-1 marks a voxel with no data** -- about half of them, because the grid box is
bigger than the brain. Those must be dropped, not read as zero.

A few experiments come back in a box of their own instead (68 x 40 x 50,
73 x 41 x 53), carrying `Offset = 0 0 0` like all the others, so the header
gives no way to place them against the atlas. Cropping them to fit would assign
every voxel to the wrong structure without complaining, so they are dropped and
listed, with the reason, in gene_drops.csv beside the table.

Orientation, checked rather than assumed. MetaImage stores x fastest, so numpy
reads the file as (z, y, x) and `transpose(2, 1, 0)` puts it in (AP, DV, ML),
the order everything else here uses. That was established two independent ways:
by matching the grid's valid-data mask against the atlas brain mask over all 48
axis permutations and flips (the winner is unambiguous, Dice 0.81 against 0.69
for the next distinct arrangement), and by the old MATLAB route, which declares
the same layout from the other side of the column-major divide. The mask test
cannot resolve a left-right flip, because the brain is nearly symmetric -- it
does not matter here, since a structure mean pools both hemispheres.

Aggregation. Voxels are assigned to structures with the CCF annotation sampled
onto the same 200 um grid, and structures are keyed by NAME, summing over the
layer-level indices, exactly as v2_region_plot does for the nano side. Two means
are computed per structure:

  full    every valid voxel inside the structure
  eroded  the same after peeling one 200 um voxel off the border

ISH is coarse and registration is imperfect, so a small structure's mean is
partly its neighbours'. The eroded version is the guard against that; both are
written so the comparison downstream can show whether it changes anything
instead of assuming it does.

Output: data\\adult_v2\\ish\\gene_region_table.csv, long format, one row per
gene and structure, plus a coverage column so a structure measured from three
voxels can be told from one measured from three thousand.

  D:\\sep_histology\\code\\tools\\venv_atlas\\Scripts\\python.exe v2_ish_regions.py [gene ...]
"""

import csv
import os
import sys
import time
from collections import defaultdict

import numpy as np
import nibabel as nib
from scipy.ndimage import binary_erosion

DATA = r'D:\sep_histology\data'
ISH_DIR = os.path.join(DATA, 'atlas_ish')
PANEL = os.path.join(DATA, 'gene_targets.csv')
CSV_MAP = os.path.join(DATA, 'atlas', 'parcellation_to_parcellation_term_membership.csv')
OUT = os.path.join(DATA, 'adult_v2', 'ish')

GRID_UM = 200
MISSING = -1.0          # the Allen flag for "no data here"
MIN_VOXELS = 3          # a structure needs this many valid 200 um voxels to get a value
GRID_DIMS = (67, 41, 58)   # the shared reference box, in the header's (x, y, z) order


class NotReferenceGrid(Exception):
    """The experiment was gridded in a box of its own, so it cannot be placed."""


def read_energy(experiment_id):
    """One Allen grid as (AP, DV, ML) at 200 um, with missing voxels as NaN.

    The header is read rather than trusted: if a future download has a different
    size or spacing, this raises instead of quietly reshaping into nonsense.
    """
    stem = os.path.join(ISH_DIR, str(experiment_id) + '_energy')
    hdr = {}
    with open(stem + '.mhd') as fh:
        for line in fh:
            if '=' in line:
                k, v = line.split('=', 1)
                hdr[k.strip()] = v.strip()
    dims = [int(x) for x in hdr['DimSize'].split()]           # (x, y, z) = (AP, DV, ML)
    spacing = {float(x) for x in hdr['ElementSpacing'].split()}
    if spacing != {float(GRID_UM)}:
        raise SystemExit(f'{experiment_id}: spacing {spacing} um, expected {GRID_UM}')
    if tuple(dims) != GRID_DIMS:
        raise NotReferenceGrid(f'grid is {tuple(dims)}, not {GRID_DIMS}')
    vol = np.fromfile(stem + '.raw', dtype=np.float32)
    if vol.size != np.prod(dims):
        raise SystemExit(f'{experiment_id}: {vol.size} values, header says {np.prod(dims)}')
    vol = vol.reshape(dims[::-1]).transpose(2, 1, 0)          # x fastest -> (AP, DV, ML)
    return np.where(vol == MISSING, np.nan, vol)


def annotation_200():
    """CCF labels on the ISH grid: the 10 um annotation sampled every 20th voxel.

    The result is one voxel smaller than the Allen grid in each axis (66 x 40 x 57
    against 67 x 41 x 58) because their box is slightly larger; the offset that
    aligns them is zero, which is what the orientation check measured.
    """
    ann = np.asarray(nib.load(os.path.join(DATA, 'atlas', 'annotation_10.nii.gz')).dataobj)
    return ann[::20, ::20, ::20]


def structure_names():
    """parcellation_index -> structure name, the same key the nano table uses."""
    names = {}
    with open(CSV_MAP, newline='', encoding='utf-8') as fh:
        for row in csv.DictReader(fh):
            if row['parcellation_term_set_name'] == 'structure':
                names[int(row['parcellation_index'])] = row['parcellation_term_name']
    return names


def region_means(vol, ann, names, eroded_ann):
    """{structure name: (full mean, eroded mean, n valid voxels, n voxels in structure)}.

    Summed over the layer-level indices that share a structure name, as on the
    nano side, so the two tables join on the same key.
    """
    out = {}
    acc = defaultdict(lambda: [0.0, 0, 0.0, 0, 0])   # sum, n, sum_ero, n_ero, n_total
    valid = np.isfinite(vol)
    for idx in np.unique(ann):
        if idx == 0 or idx not in names:
            continue
        m = ann == idx
        mv = m & valid
        a = acc[names[idx]]
        a[0] += float(vol[mv].sum()); a[1] += int(mv.sum()); a[4] += int(m.sum())
        me = (eroded_ann == idx) & valid
        a[2] += float(vol[me].sum()); a[3] += int(me.sum())
    for name, (s, n, se, ne, ntot) in acc.items():
        if n >= MIN_VOXELS:
            out[name] = (s / n, (se / ne) if ne >= MIN_VOXELS else np.nan, n, ntot)
    return out


def main(only=None):
    os.makedirs(OUT, exist_ok=True)
    panel = [r for r in csv.DictReader(open(PANEL, newline='', encoding='utf-8'))]
    if only:
        want = {g.lower() for g in only}
        panel = [r for r in panel if r['symbol'].lower() in want]
    names = structure_names()
    ann_full = annotation_200()

    # erode each structure by one voxel, once, and reuse it for every gene
    print('eroding structure masks once (200 um, one voxel)...', flush=True)
    eroded = np.zeros_like(ann_full)
    for idx in np.unique(ann_full):
        if idx == 0:
            continue
        m = ann_full == idx
        e = binary_erosion(m)
        eroded[e if e.any() else m] = idx        # keep structures too small to erode

    rows, dropped = [], []
    for i, gene in enumerate(panel, 1):
        t0 = time.time()
        sym, eid = gene['symbol'], gene['experiment_id']
        try:
            vol = read_energy(eid)
        except FileNotFoundError:
            dropped.append(dict(symbol=sym, experiment_id=eid, reason='grid not downloaded'))
            print(f'{i:3d}/{len(panel)} {sym:10s} DROPPED  grid not downloaded', flush=True)
            continue
        except NotReferenceGrid as why:
            dropped.append(dict(symbol=sym, experiment_id=eid, reason=str(why)))
            print(f'{i:3d}/{len(panel)} {sym:10s} DROPPED  {why}', flush=True)
            continue
        vol = vol[:ann_full.shape[0], :ann_full.shape[1], :ann_full.shape[2]]
        rm = region_means(vol, ann_full, names, eroded)
        for name, (full, ero, n, ntot) in rm.items():
            rows.append(dict(symbol=sym, experiment_id=eid, category=gene['category'],
                             structure=name, ish_mean=f'{full:.6g}',
                             ish_mean_eroded=('' if np.isnan(ero) else f'{ero:.6g}'),
                             n_voxels=n, n_voxels_structure=ntot,
                             coverage=f'{n / max(ntot, 1):.3f}'))
        print(f'{i:3d}/{len(panel)} {sym:10s} {len(rm):3d} structures   {time.time() - t0:.1f} s', flush=True)

    path = os.path.join(OUT, 'gene_region_table.csv')
    with open(path, 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(f'\n{len(rows):,} rows, {len(panel) - len(dropped)} genes -> {path}',
          flush=True)

    path = os.path.join(OUT, 'gene_drops.csv')
    with open(path, 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=['symbol', 'experiment_id', 'reason'])
        w.writeheader(); w.writerows(dropped)
    print(f'{len(dropped)} genes dropped -> {path}', flush=True)


if __name__ == '__main__':
    main(sys.argv[1:] or None)
