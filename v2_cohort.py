"""
v2, step 3: cohort volumes, built in the adult CCF from the per-mouse files.

Everything is averaged on one grid, 660 x 400 x 570 at 20 um, because v2_to_ccf
has already carried each brain there: the young through the DeMBA -> CCF
transform of its own age, the adults by placement alone. That is what lets a
pooled young group mix P20 and P16 without either age being carried by the
other's deformation field.

Five readings of the same background-subtracted signal, the same five the
region tables carry, so a number in a table and a colour in a map mean the
same thing. Two divide by a CHANNEL, voxel by voxel; three divide by a NUMBER
measured on the brain itself:
  ratio   sig / auto per voxel (auto smoothed by one 20 um voxel so a dark
          voxel cannot blow it up): nano per unit autofluorescence, the
          internal standard.
          Read this one carefully. Measured in the isocortex, the young brains
          sit 2.0 log2 below the adults in nano and 1.0 log2 below them in
          auto, so the 1.0 log2 that survives in the ratio is what is left
          after dividing one age-dependent quantity by another. It is not an
          absolute measurement: autofluorescence rises with age (lipofuscin,
          tissue density), so this reading understates a real pup deficit and
          would overstate a pup excess. Treat it as a bound, not a value.
  sepratio
          sig / SEP per voxel, the same arithmetic with the green channel as
          the denominator. SEP is the tag on GluA1 itself, so this one reads as
          receptor on the membrane per unit receptor expressed -- a quantity
          with a meaning, where nano/auto is only nano per unit tissue. It has
          its own weakness, though, and the opposite one: if expression itself
          changes with age, this reading divides the effect out. Between the
          two, ratio is blind to expression and sepratio is blind to a common
          scaling of both pools; they are kept side by side for that reason.
  cref    sig / (that mouse's isocortex mean of sig, measured natively before
          any warp): a pure scale, so region ratios within a mouse survive
          exactly. Cortex-relative by construction, hence blind to a change
          that moves the whole cortex.
  subref  sig / (that mouse's subcortex mean, excluding HPF and STR): the same
          idea with a reference that is neither the cortex nor the two
          structures that dominate the scale.
  zref    range-matched: log2(sig / cortex mean), minus that brain's median
          over structures and divided by its own p90-p10 spread. Level AND
          dynamic range are then the same in every brain, which matters
          because the pup brain is genuinely flatter (spread 0.89 +- 0.25 log2
          against 1.79 +- 0.24). It is the only reading here that is not
          linear in the signal, so its maps are a position within a range, not
          an intensity; and the compression it removes may itself be the
          finding, so it is read beside cref, not instead of it.

Per voxel the cohort gets the mean over the mice that have tissue there, the
SD and the count. No voxel is required to be covered by every mouse; the n map
says what each mean rests on and the reader thresholds it.

Cohorts: young (P16 + P20 + P22 pooled, the group Sami wants), young_P20 (the
sensitivity check with the P20 brains alone), young_P16, young_P22, naive, rws,
and adult = naive + rws.

Output: data/comparisons_v2/ccf/<cohort>/{ratio,cref}_{mean,sd,n}.npy + mice.txt

  D:\\sep_histology\\code\\tools\\venv_atlas\\Scripts\\python.exe v2_cohort.py
"""

import os

import numpy as np
from scipy.ndimage import gaussian_filter

from v2_per_mouse import DATA, annotation_20, MICE, CSV_MAP, OUT as PER_MOUSE

V2 = os.path.join(DATA, 'comparisons_v2')
PER_MOUSE_CCF = os.path.join(V2, 'per_mouse_ccf')
OUT_ROOT = os.path.join(V2, 'ccf')
RATIO_CLIP = 20.0
Z_FLOOR = 0.02        # of the cortex mean, so log2 stays finite in the dimmest tissue
MODES = ('ratio', 'sepratio', 'cref', 'subref', 'zref')
NOT_SUBCORTEX = {'Isocortex', 'HPF', 'STR', 'OLF', 'CTXsp', 'fiber tracts', 'VS', 'CB', ''}

YOUNG_P20 = ['MG897_SepGluA_P20', 'MG903_SepGluA_P20', 'MG913_SepGluA_P20',
             'MG909_SepGluA_P20', 'MG910_SepGluA_P20']
YOUNG_P16 = ['MG911_SepGluA_P16']
YOUNG_P22 = ['MG904_SepGluA_P22']
NAIVE = ['CGF027_Gria1', 'CGF028_Gria1', 'CGF033_Gria1', 'CGF034_Gria1', 'CGF035_Gria1']
RWS = ['MG691_Gria1', 'MG692_Gria1', 'MG693_Gria1', 'MG736_Gria1', 'MG737_Gria1']
COHORTS = {
    'young': YOUNG_P20 + YOUNG_P16 + YOUNG_P22,
    'young_P20': YOUNG_P20,
    'young_P16': YOUNG_P16,
    'young_P22': YOUNG_P22,
    'naive': NAIVE,
    'rws': RWS,
    'adult': NAIVE + RWS,
}


def mouse_scalars(mouse):
    """The per-brain numbers every reading divides by, measured on the brain's OWN
    atlas before any warp, and cached beside the per-mouse file.

    cortex_mean and subcortex_mean are plain means over their labels; z_median
    and z_spread describe that brain's distribution ACROSS STRUCTURES of
    log2(structure mean / cortex mean), which is what the range match needs.
    Structure level, not voxel level, so a large structure cannot set the
    spread on its own.
    """
    src = os.path.join(PER_MOUSE, mouse + '.npz')
    cache = os.path.join(PER_MOUSE, mouse + '_scalars.npz')
    # The cache carries the modification time of the file it was computed from,
    # so re-running v2_per_mouse silently invalidates it. Without that a changed
    # tissue mask would go on being divided by the old cortex mean, and nothing
    # on screen would say so.
    if os.path.exists(cache):
        z = np.load(cache)
        if 'src_mtime' in z.files and float(z['src_mtime']) == os.path.getmtime(src):
            return {k: float(z[k]) for k in z.files if k != 'src_mtime'}
    import csv as _csv
    stru, divi = {}, {}
    with open(CSV_MAP, newline='', encoding='utf-8') as fh:
        for r in _csv.DictReader(fh):
            i = int(r['parcellation_index'])
            if r['parcellation_term_set_name'] == 'structure':
                stru[i] = r['parcellation_term_acronym']
            elif r['parcellation_term_set_name'] == 'division':
                divi[i] = r['parcellation_term_acronym']
    ann = annotation_20(MICE[mouse][1])
    z = np.load(os.path.join(PER_MOUSE, mouse + '.npz'))
    sig = z['sig'].astype(np.float32); tissue = z['tissue']
    lab = ann[tissue]; nlab = int(ann.max()) + 1
    n = np.bincount(lab, minlength=nlab)
    tot = np.bincount(lab, weights=sig[tissue], minlength=nlab)
    iso = [i for i in stru if divi.get(i) == 'Isocortex' and i < nlab]
    sub = [i for i in stru if divi.get(i, '') not in NOT_SUBCORTEX and i < nlab]
    cortex_mean = tot[iso].sum() / n[iso].sum()
    sub_mean = tot[sub].sum() / n[sub].sum()
    per_struct = {}
    for i in np.nonzero(n)[0]:
        if i == 0 or i not in stru:
            continue
        a, b = per_struct.get(stru[i], (0, 0.0))
        per_struct[stru[i]] = (a + int(n[i]), b + tot[i])
    vals = np.array([np.log2(t / c / cortex_mean) for c, t in per_struct.values()
                     if c >= 250 and t > 0])
    p10, med, p90 = np.percentile(vals, [10, 50, 90])
    out = dict(cortex_mean=float(cortex_mean), subcortex_mean=float(sub_mean),
               z_median=float(med), z_spread=float(max(p90 - p10, 1e-6)))
    np.savez(cache, src_mtime=os.path.getmtime(src), **out)
    return out


def mouse_modes(mouse):
    """(dict of mode -> volume with NaN off tissue, tissue mask) for one brain in CCF."""
    z = np.load(os.path.join(PER_MOUSE_CCF, mouse + '.npz'))
    sig = z['sig'].astype(np.float32); auto = z['auto'].astype(np.float32); tissue = z['tissue']
    if 'sep' not in z.files:
        raise SystemExit(f'{mouse}: no SEP channel in its per-mouse CCF file. Run\n'
                         f'  P4bis_add_sep_channel.m for this brain, then v2_per_mouse.py and v2_to_ccf.py.')
    sep = z['sep'].astype(np.float32)
    sc = mouse_scalars(mouse)

    def per_unit(ref):
        """sig divided by a reference CHANNEL, voxel by voxel.

        The denominator is smoothed by one 20 um voxel first, so a single dark
        voxel in the reference cannot blow the ratio up; the smoothing is
        normalised by the mask, so tissue at the edge is not divided by the
        black outside it.
        """
        ref_s = gaussian_filter(np.where(tissue, ref, 0), 1.0) / np.maximum(
            gaussian_filter(tissue.astype(np.float32), 1.0), 1e-3)
        r = np.where(tissue & (ref_s > 0), sig / np.maximum(ref_s, 1e-3), np.nan)
        return np.clip(r, -RATIO_CLIP, RATIO_CLIP).astype(np.float32)

    cref = np.where(tissue, sig / sc['cortex_mean'], np.nan)
    subref = np.where(tissue, sig / sc['subcortex_mean'], np.nan)
    floored = np.maximum(sig / sc['cortex_mean'], Z_FLOOR)
    zref = np.where(tissue, (np.log2(floored) - sc['z_median']) / sc['z_spread'], np.nan)
    return ({'ratio': per_unit(auto), 'sepratio': per_unit(sep),
             'cref': cref.astype(np.float32), 'subref': subref.astype(np.float32),
             'zref': zref.astype(np.float32)}, tissue)


def main():
    for cohort, mice in COHORTS.items():
        out = os.path.join(OUT_ROOT, cohort); os.makedirs(out, exist_ok=True)
        acc = None
        for mouse in mice:
            modes, tissue = mouse_modes(mouse)
            if acc is None:
                acc = {k: [np.zeros(v.shape, np.float64), np.zeros(v.shape, np.float64), np.zeros(v.shape, np.int16)]
                       for k, v in modes.items()}
            for k, v in modes.items():
                w = np.nan_to_num(v)
                acc[k][0] += w; acc[k][1] += w * w; acc[k][2] += tissue
        for k, (s, ss, n) in acc.items():
            nf = np.maximum(n, 1).astype(np.float64)
            mean = np.where(n > 0, s / nf, np.nan).astype(np.float32)
            var = np.where(n > 1, (ss - s * s / nf) / np.maximum(nf - 1, 1), np.nan)
            np.save(os.path.join(out, f'{k}_mean.npy'), mean)
            np.save(os.path.join(out, f'{k}_sd.npy'), np.sqrt(np.maximum(var, 0)).astype(np.float32))
            np.save(os.path.join(out, f'{k}_n.npy'), n)
        with open(os.path.join(out, 'mice.txt'), 'w') as fh:
            fh.write('\n'.join(mice) + '\n')
        n = acc['cref'][2]
        print(f'{cohort:10s} {len(mice):2d} mice | voxels with n>=1 {int((n > 0).sum()):>11,d} | '
              f'with all mice {int((n == len(mice)).sum()):>11,d}', flush=True)


if __name__ == '__main__':
    main()
