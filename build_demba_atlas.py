"""
Build the DeMBA atlas folder for one postnatal age, the way atlas_demba_p20 was
built, so a brain of that age can be registered to a template of its own age.

BrainGlobe publishes DeMBA (Carey 2025) with Allen CCFv3 segmentation at every
age from P4 to P56, 20 um isotropic, oriented (AP, DV, ML) like the volumes this
pipeline indexes. Three things still have to happen before MATLAB can use one:

  names      LightSuite finds the atlas with which('average_template_10.nii.gz')
             in fourteen vendored files, so the files must carry the _10 name
             even though they hold 20 um data. The real resolution lives in
             get_atlas (res_um) and in each mouse's local_settings.txt.
  ID space   BrainGlobe ships Allen *structure IDs*; everything downstream here
             (get_allen_region_mask, P5, P8) resolves regions through
             *parcellation_index*. The two collide numerically without meaning
             the same thing, so the annotation is remapped. The original is kept
             beside it as annotation_structureids_original.nii.gz.
  AP crop    the adult brains are cropped to CCF planes [180 1079]; the same
             anatomy on a developmental grid is not that range divided by two,
             because DeMBA is about 11% longer in AP for the same brain. The
             crop is measured here by two independent methods (see below) and
             written to aplims.txt, which get_atlas reads.

The two crop measurements are the same ones used for P20, where they agreed
([63 559] and [62 562]):

  area profile   the brain's cross-sectional area along AP has a characteristic
                 shape (bulb, the cortex/hippocampus bulge, the taper into
                 brainstem). The crop whose normalised profile best matches the
                 adult crop's selects the same anatomy. Not circular: it never
                 looks at image intensities.
  region COM     regress the AP centre of mass of every region present in both
                 atlases, CCF_plane = m * DeMBA_plane + b, then invert it at the
                 adult crop limits. Independent of the shape of the brain as a
                 whole, and its slope also measures the AP stretch.

  D:\\sep_histology\\code\\tools\\venv_atlas\\Scripts\\python.exe build_demba_atlas.py 16
  ... then add the printed case to get_atlas.m if it is not there yet.
"""

import csv
import re
import sys
from pathlib import Path

import nibabel as nib
import numpy as np

DATA = Path(r'D:\sep_histology\data')
CCF_ANN = DATA / 'atlas' / 'annotation_10.nii.gz'
PARCELLATION = DATA / 'atlas' / 'parcellation.csv'
CCF_LIMS = (180, 1079)          # the adult crop, 1-based inclusive, 10 um planes
MIN_VOX_CCF, MIN_VOX_DEMBA = 5000, 600      # a region has to be solid in both to be fitted


def area_profile(ann):
    """Cross-sectional area of labelled brain per AP plane."""
    return (ann > 0).sum(axis=(1, 2)).astype(float)


def crop_by_area_profile(dem_ann, ccf_ann):
    """The DeMBA crop whose normalised area profile best matches the adult crop's."""
    ccf_area = area_profile(ccf_ann)[CCF_LIMS[0] - 1:CCF_LIMS[1]]
    ref = ccf_area / ccf_area.max()
    x = np.linspace(0, 1, 400)
    ref_i = np.interp(x, np.linspace(0, 1, len(ref)), ref)

    dem_area = area_profile(dem_ann)
    nz = np.flatnonzero(dem_area)
    d_first, d_last = int(nz[0]), int(nz[-1])

    def rms(lo, hi):
        seg = dem_area[lo - 1:hi]
        if seg.size < 50 or seg.max() == 0:
            return np.inf
        seg_i = np.interp(x, np.linspace(0, 1, seg.size), seg / seg.max())
        return float(np.sqrt(np.mean((ref_i - seg_i) ** 2)))

    best, step = None, 4
    los = range(max(1, d_first - 20), d_first + 140, step)
    his = range(d_last - 140, min(len(dem_area), d_last + 20), step)
    for lo in los:
        for hi in his:
            if hi - lo < 100:
                continue
            r = rms(lo, hi)
            if best is None or r < best[0]:
                best = (r, lo, hi)
    for _ in range(2):                       # refine around the coarse optimum
        step = max(1, step // 2)
        r0, lo0, hi0 = best
        for lo in range(lo0 - 2 * step, lo0 + 2 * step + 1, step):
            for hi in range(hi0 - 2 * step, hi0 + 2 * step + 1, step):
                if hi - lo < 100:
                    continue
                r = rms(lo, hi)
                if r < best[0]:
                    best = (r, lo, hi)
    return best[1], best[2], best[0], (d_first + 1, d_last + 1)


def crop_by_region_com(dem_idx, ccf_ann):
    """Invert the AP regression over regions shared by the two atlases."""
    shared = np.intersect1d(np.unique(ccf_ann), np.unique(dem_idx))
    shared = shared[shared != 0]

    def com(vol):
        n = len(shared)
        w, c = np.zeros(n), np.zeros(n)
        for i in range(vol.shape[0]):
            vals, cnts = np.unique(vol[i], return_counts=True)
            p = np.searchsorted(shared, vals)
            p[p >= n] = 0
            ok = shared[p] == vals
            np.add.at(c, p[ok], cnts[ok])
            np.add.at(w, p[ok], cnts[ok] * i)
        with np.errstate(invalid='ignore', divide='ignore'):
            return w / c, c

    com_c, tot_c = com(ccf_ann)
    com_d, tot_d = com(dem_idx)
    keep = (tot_c >= MIN_VOX_CCF) & (tot_d >= MIN_VOX_DEMBA)
    m, b = np.polyfit(com_d[keep], com_c[keep], 1)
    resid = com_c[keep] - (m * com_d[keep] + b)
    lo, hi = round((CCF_LIMS[0] - b) / m), round((CCF_LIMS[1] - b) / m)
    return int(lo), int(hi), float(m), float(resid.std()), int(keep.sum())


def main(age):
    name = f'demba_allen_seg_dev_mouse_p{age}_20um'
    out_dir = DATA / f'atlas_demba_p{age}'
    print(f'fetching {name} (downloads on first use) ...', flush=True)
    from brainglobe_atlasapi import BrainGlobeAtlas
    atlas = BrainGlobeAtlas(name)
    tmpl = np.asarray(atlas.template)
    ann = np.asarray(atlas.annotation)
    res = atlas.metadata['resolution']
    assert tuple(res) == (20, 20, 20), f'expected 20 um isotropic, got {res}'
    assert atlas.orientation == 'asr', f'expected asr (AP, DV, ML), got {atlas.orientation}'
    print(f'  shape {ann.shape}, orientation {atlas.orientation}, {len(np.unique(ann))} labels', flush=True)

    # ---------------------------------------------- structure ids -> parcellation_index
    sid_to_index = {}
    with open(PARCELLATION, newline='', encoding='utf-8') as fh:
        for row in csv.DictReader(fh):
            m = re.search(r'AllenCCF-Annotation-\d+-(\d+)$', row['label'])
            if m:
                sid_to_index[int(m.group(1))] = int(row['parcellation_index'])
    present = np.unique(ann)
    lut = np.array([0 if v == 0 else sid_to_index.get(int(v), 0) for v in present], dtype=np.uint16)
    missing = [int(v) for v, t in zip(present, lut) if v != 0 and t == 0]
    ann_idx = lut[np.searchsorted(present, ann)].astype(np.uint16)
    kept = (ann_idx != 0).sum() / max((ann != 0).sum(), 1)
    print(f'  ids {len(present)}, untranslatable {len(missing)}, labelled voxels kept {kept:.4f}', flush=True)
    if kept < 0.999:
        raise SystemExit(f'refusing to write: {1 - kept:.3%} of labelled voxels lost in the id remap')

    # ------------------------------------------------------------------- crop
    ccf_ann = np.asarray(nib.load(str(CCF_ANN)).dataobj)
    lo_a, hi_a, rms, brain = crop_by_area_profile(ann_idx, ccf_ann)
    print(f'  brain spans AP planes {brain[0]}..{brain[1]} of {ann.shape[0]}', flush=True)
    print(f'  crop, area profile : [{lo_a} {hi_a}]  ({hi_a - lo_a + 1} planes, '
          f'{(hi_a - lo_a + 1) * 20 / 1000:.2f} mm, profile RMS {rms:.4f})', flush=True)
    lo_r, hi_r, slope, sd, n_reg = crop_by_region_com(ann_idx, ccf_ann)
    print(f'  crop, region COM   : [{lo_r} {hi_r}]  (from {n_reg} regions, slope {slope:.3f} '
          f'CCF planes per DeMBA plane, residual {sd:.1f} CCF planes = {sd * 10 / 1000:.3f} mm)', flush=True)
    if abs(lo_a - lo_r) > 25 or abs(hi_a - hi_r) > 25:
        print('  WARNING: the two methods disagree by more than 25 planes (0.5 mm). '
              'Check before registering anything to this atlas.', flush=True)

    # ----------------------------------------------------------------- write
    out_dir.mkdir(parents=True, exist_ok=True)
    affine = np.diag([0.02, 0.02, 0.02, 1.0])      # 20 um isotropic, mm units
    for arr, fname, dtype in ((tmpl, 'average_template_10.nii.gz', np.uint16),
                              (ann, 'annotation_structureids_original.nii.gz', np.uint32),
                              (ann_idx, 'annotation_10.nii.gz', np.uint16)):
        img = nib.Nifti1Image(arr.astype(dtype), affine)
        img.header.set_xyzt_units('mm')
        nib.save(img, str(out_dir / fname))
        print(f'  wrote {out_dir / fname}  ({(out_dir / fname).stat().st_size / 1e6:.1f} MB)', flush=True)
    (out_dir / 'aplims.txt').write_text(f'{lo_a} {hi_a}\n')
    (out_dir / 'source.txt').write_text(
        f'DeMBA P{age}, Allen CCFv3 segmentation, 20 um isotropic.\n'
        f'Downloaded via brainglobe-atlasapi as {name}.\n'
        'Carey et al. 2025, Nat Commun 16:8108. doi:10.1101/2024.06.14.598876\n'
        f'grid shape (AP, DV, ML): {ann.shape}, orientation asr\n'
        f'brain spans AP planes {brain[0]}..{brain[1]}\n'
        f'annotation remapped to Allen parcellation_index ({len(present)} ids, '
        f'{len(missing)} untranslatable, {kept:.4f} of labelled voxels kept)\n'
        f'crop equivalent to the adult CCF [180 1079]:\n'
        f'  area profile : [{lo_a} {hi_a}]  (RMS {rms:.4f})\n'
        f'  region COM   : [{lo_r} {hi_r}]  ({n_reg} regions, slope {slope:.3f}, '
        f'residual {sd:.1f} CCF planes)\n'
        f'aplims.txt holds the area-profile value, which is what get_atlas uses.\n'
        f'Built by build_demba_atlas.py.\n')
    print(f'\nwrote {out_dir / "aplims.txt"} ({lo_a} {hi_a}) and source.txt')
    print(f"get_atlas('demba_p{age}') will pick this up; P4 needs atlas_key = 'demba_p{age}'.")


if __name__ == '__main__':
    if len(sys.argv) != 2 or not sys.argv[1].isdigit():
        raise SystemExit(__doc__.strip().splitlines()[-2].strip())
    main(int(sys.argv[1]))
