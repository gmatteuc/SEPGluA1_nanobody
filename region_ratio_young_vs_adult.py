"""
P20 vs adult, region by region, under several choices of reference and of
what "signal" means.

Reads the per-mouse raw tables written by region_means_raw_per_mouse.py (nano,
and auto if present) and turns every region into log2(region / reference) PER
MOUSE, so the cohorts can be compared with a spread and a test (P20 n = 3,
adults n = 10).

Three definitions of the signal, each a MODE:
  nano          the raw nano mean, as is
  nano-bg       raw nano minus the mouse's off-tissue background. An additive
                offset compresses every ratio toward 1, and it compresses the
                pups' more because their raw levels are 2-3x lower, so this is
                the one to trust for ratios.
  nano/auto     (nano - bg_nano) / (auto - bg_auto): the autofluorescence
                channel as an internal standard, per voxel. The only mode with
                a 'none' reference, i.e. an absolute-ish comparison, valid to
                the extent that autofluorescence itself does not change with
                age (it is the assumption, not a fact).

And the REFERENCES. Each answers a different question and none settles the
critical-period prediction on its own:
  isocortex        is VIS/SS a larger share OF THE CORTEX in pups? Blind to a
                   whole-cortex upregulation.
  thalamus / VPM+VPL
                   is cortex up relative to its input nuclei?
  pallidum, hypothalamus, midbrain
                   candidate "stable" subcortical anchors. None is known to be
                   developmentally invariant for surface GluA1; the spread
                   between them is the honest uncertainty of the answer.
  subcortex-HPF-STR
                   everything below the cortex except the two structures that
                   dominate the scale at both ages.
  brain            all tissue.
  none             the value itself (nano/auto mode only).

Alongside P20-vs-adult, naive-vs-rws is printed under the same reference: two
adult cohorts with no developmental difference between them. A P20 effect
that is not clearly larger than that one is not an effect.

Output: region_ratio_young_vs_adult.csv (all regions x modes x references),
a per-mouse reference table, and cortex tables per mode.

  D:\\sep_histology\\code\\tools\\venv_atlas\\Scripts\\python.exe region_ratio_young_vs_adult.py
"""

import csv
import math
import os
from collections import defaultdict

import numpy as np

DATA = r'D:\sep_histology\data'
OUT = os.path.join(DATA, 'comparisons', 'young_P20_vs_adult_nano')
TABLES = {'nano': os.path.join(OUT, 'region_means_raw_per_mouse.csv'),
          'auto': os.path.join(OUT, 'region_means_raw_per_mouse_auto.csv')}
MIN_VOX = 2000

NOT_SUBCORTEX = {'Isocortex', 'HPF', 'STR', 'OLF', 'CTXsp', 'fiber tracts', 'VS', 'CB', '', 'bg'}
REFERENCES = {
    'isocortex':          lambda n, a, d: d == 'Isocortex',
    'thalamus':           lambda n, a, d: d == 'TH',
    'VPM+VPL':            lambda n, a, d: a in ('VPM', 'VPL'),
    'pallidum':           lambda n, a, d: d == 'PAL',
    'hypothalamus':       lambda n, a, d: d == 'HY',
    'midbrain':           lambda n, a, d: d == 'MB',
    'subcortex-HPF-STR':  lambda n, a, d: d not in NOT_SUBCORTEX,
    'brain':              lambda n, a, d: d not in ('fiber tracts', 'VS', '', 'bg'),
    'none':               None,
}
CORTEX_SHOW = ('VISp', 'VISl', 'VISal', 'VISrl', 'VISpm', 'VISam', 'VISli',
               'SSp-bfd', 'SSp-ll', 'SSp-ul', 'SSp-m', 'SSp-n', 'SSp-tr', 'SSp-un', 'SSs',
               'AUDp', 'AUDd', 'AUDv', 'MOp', 'MOs', 'RSPagl', 'RSPd', 'RSPv', 'ACAd', 'ACAv',
               'ORBl', 'ORBm', 'ORBvl', 'PL', 'ILA', 'TEa', 'ECT', 'PERI', 'AId', 'AIv', 'AIp', 'GU', 'VISC')


def welch(a, b):
    """Welch t and two-sided p for two samples."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 2 or len(b) < 2:
        return float('nan'), float('nan')
    va, vb = a.var(ddof=1) / len(a), b.var(ddof=1) / len(b)
    if va + vb == 0:
        return float('nan'), float('nan')
    t = (a.mean() - b.mean()) / math.sqrt(va + vb)
    df = (va + vb) ** 2 / (va ** 2 / (len(a) - 1) + vb ** 2 / (len(b) - 1))
    try:
        from scipy.stats import t as tdist
        p = 2 * tdist.sf(abs(t), df)
    except ImportError:
        p = float('nan')
    return t, p


def load(path):
    """-> per_mouse[mouse][structure] = (sum, count), bg[mouse], cohort_of, meta[structure] = (acronym, division)"""
    per_mouse = defaultdict(lambda: defaultdict(lambda: [0.0, 0]))
    bg, cohort_of, meta = {}, {}, {}
    with open(path, newline='', encoding='utf-8') as fh:
        for r in csv.DictReader(fh):
            m = r['mouse']; cohort_of[m] = r['cohort']
            if r['acronym'] == 'bg':
                bg[m] = float(r['mean_raw']); continue
            n = r['structure']; meta[n] = (r['acronym'], r['division'])
            cell = per_mouse[m][n]
            cell[0] += float(r['mean_raw']) * int(r['n_vox']); cell[1] += int(r['n_vox'])
    return per_mouse, bg, cohort_of, meta


def main():
    nano, bg_n, cohort_of, meta = load(TABLES['nano'])
    have_auto = os.path.exists(TABLES['auto'])
    if have_auto:
        auto, bg_a, _, _ = load(TABLES['auto'])
    mice = sorted(nano, key=lambda m: (cohort_of[m], m))
    cohorts = {c: [m for m in mice if cohort_of[m] == c] for c in ('young_P20', 'naive', 'rws')}
    adults = cohorts['naive'] + cohorts['rws']

    # value[mode][mouse][structure] -> per-voxel-mean signal, and its voxel count
    def region_value(mode, m, n):
        sv, cv = nano[m].get(n, (0.0, 0))
        if cv < MIN_VOX:
            return None, 0
        v = sv / cv
        if mode == 'nano':
            return v, cv
        v -= bg_n[m]
        if mode == 'nano-bg':
            return v, cv
        sa, ca = auto[m].get(n, (0.0, 0))
        if ca < MIN_VOX:
            return None, 0
        a = sa / ca - bg_a[m]
        return (v / a if a > 0 else None), cv

    modes = ['nano', 'nano-bg'] + (['nano/auto'] if have_auto else [])
    structures = sorted(meta, key=lambda n: (meta[n][1], meta[n][0]))
    value = {mode: {m: {n: region_value(mode, m, n) for n in structures} for m in mice} for mode in modes}

    # reference level per mouse and mode: voxel-weighted mean of the member structures
    ref_level = {}
    for mode in modes:
        for ref, pred in REFERENCES.items():
            for m in mice:
                if pred is None:
                    ref_level[(mode, ref, m)] = 1.0; continue
                s = c = 0.0
                for n in structures:
                    v, cv = value[mode][m][n]
                    if v is not None and pred(n, *meta[n]):
                        s += v * cv; c += cv
                ref_level[(mode, ref, m)] = s / c if c else float('nan')

    print('per mouse: background (off-tissue median) and reference levels, nano-bg mode (raw units, comparable only within a cohort)')
    refs_show = [r for r in REFERENCES if r != 'none']
    print(f'  {"mouse":20s} {"cohort":10s} {"bg nano":>8s} ' + (f'{"bg auto":>8s} ' if have_auto else '') + ' '.join(f'{r:>17s}' for r in refs_show))
    for m in mice:
        print(f'  {m:20s} {cohort_of[m]:10s} {bg_n[m]:8.0f} ' + (f'{bg_a[m]:8.0f} ' if have_auto else '')
              + ' '.join(f'{ref_level[("nano-bg", r, m)]:17.0f}' for r in refs_show))
    if have_auto:
        print('  nano/auto of the isocortex per mouse (bg-subtracted both):  '
              + '  '.join(f'{m.split("_")[0]} {ref_level[("nano/auto", "isocortex", m)]:.2f}' for m in mice))

    rows = []
    for mode in modes:
        for n in structures:
            a, d = meta[n]
            for ref in REFERENCES:
                if ref == 'none' and mode != 'nano/auto':
                    continue
                v = {}
                for m in mice:
                    x, cv = value[mode][m][n]
                    rl = ref_level[(mode, ref, m)]
                    if x is not None and x > 0 and rl > 0:
                        v[m] = math.log2(x / rl)
                yp = [v[m] for m in cohorts['young_P20'] if m in v]
                ad = [v[m] for m in adults if m in v]
                nv = [v[m] for m in cohorts['naive'] if m in v]
                rw = [v[m] for m in cohorts['rws'] if m in v]
                if len(yp) < 2 or len(ad) < 4:
                    continue
                t, p = welch(yp, ad)
                t2, p2 = welch(nv, rw)
                rows.append(dict(mode=mode, reference=ref, structure=n, acronym=a, division=d,
                                 n_P20=len(yp), n_adult=len(ad),
                                 P20_mean=np.mean(yp), P20_sem=np.std(yp, ddof=1) / math.sqrt(len(yp)),
                                 adult_mean=np.mean(ad), adult_sem=np.std(ad, ddof=1) / math.sqrt(len(ad)),
                                 diff_log2=np.mean(yp) - np.mean(ad), t=t, p=p,
                                 naive_vs_rws_log2=(np.mean(nv) - np.mean(rw)) if nv and rw else float('nan'),
                                 naive_vs_rws_p=p2))
    out = os.path.join(OUT, 'region_ratio_young_vs_adult.csv')
    with open(out, 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow({k: (f'{v:.4f}' if isinstance(v, float) else v) for k, v in r.items()})
    print(f'\nwrote {out}: {len(rows)} rows')

    by = {(r['mode'], r['acronym'], r['reference']): r for r in rows}
    for mode in modes:
        refs = [r for r in REFERENCES if r != 'none' or mode == 'nano/auto']
        print(f'\n=== MODE {mode} ===  CORTEX: log2(P20 / adult) of the region relative to each reference '
              '(* p<0.05, ** p<0.01 Welch, n=3 vs 10). Last column: naive vs rws under the isocortex reference.')
        print(f'  {"area":9s} ' + ' '.join(f'{r:>17s}' for r in refs) + f' {"naive-rws":>10s}')
        for a in CORTEX_SHOW:
            cells = []
            for ref in refs:
                r = by.get((mode, a, ref))
                if r is None:
                    cells.append(f'{"--":>17s}'); continue
                star = '**' if r['p'] < 0.01 else ('*' if r['p'] < 0.05 else '')
                cells.append(f'{r["diff_log2"]:+7.2f}{star:2s}' + ' ' * 8)
            r = by.get((mode, a, 'isocortex'))
            tail = f'{r["naive_vs_rws_log2"]:+8.2f}' if r else ''
            print(f'  {a:9s} ' + ' '.join(cells) + f' {tail}')

    # brain-wide, for the trusted mode(s): what moves, with the naive-rws null next to it
    for mode, ref in (('nano-bg', 'isocortex'), ('nano-bg', 'subcortex-HPF-STR')) + ((('nano/auto', 'none'),) if have_auto else ()):
        sub = sorted([r for r in rows if r['mode'] == mode and r['reference'] == ref and r['p'] < 0.05],
                     key=lambda r: -r['diff_log2'])
        print(f'\n=== MODE {mode}, reference {ref}: significant (p<0.05) structures, MORE in P20 first ===')
        for r in sub[:15] + ([('...',)] if len(sub) > 30 else []) + sub[-15:]:
            if isinstance(r, tuple):
                print('  ...'); continue
            print(f'  {r["acronym"]:9s} {r["structure"][:36]:36s} {r["diff_log2"]:+6.2f}  p={r["p"]:.3f}  '
                  f'naive-rws {r["naive_vs_rws_log2"]:+5.2f}  [{r["division"]}]')


if __name__ == '__main__':
    main()
