"""
Download the Allen expression grids the panel asks for, once.

One zip per experiment from `api.brain-map.org/grid_data/download/<id>`, holding
`energy.mhd` and `energy.raw`, unpacked next to the ones already there in
data\\atlas_ish under the experiment id. Grids already on disk are left alone, so
this can be re-run at any time and will only fetch what is missing -- including
after a partial run, which matters because there are several hundred of them.

Two things are recorded rather than swallowed. A grid that 404s (some
experiments have no downloadable grid at all, which is not a local problem) and
a grid whose header is not the shared 67 x 41 x 58 reference box (a few
experiments come in a box of their own and cannot be placed against the atlas;
v2_ish_regions drops them for the same reason). Both go to fetch_failures.csv
with the reason, so the panel's real size is a number on disk and not a guess.

  D:\\sep_histology\\code\\tools\\venv_atlas\\Scripts\\python.exe mapping\\v2_panel_fetch.py
"""

import csv
import io
import os
import time
import urllib.error
import urllib.request
import zipfile

from v2_paths import DATA
PANEL = os.path.join(DATA, 'adult_v2', 'panel', 'panel_v2.csv')
DEST = os.path.join(DATA, 'atlas_ish')
OUT = os.path.join(DATA, 'adult_v2', 'panel')

URL = 'http://api.brain-map.org/grid_data/download/{}'
GRID_DIMS = (67, 41, 58)
TIMEOUT = 180
RETRIES = 2


def already_there(eid):
    return all(os.path.exists(os.path.join(DEST, f'{eid}_energy{ext}'))
               for ext in ('.mhd', '.raw'))


def dims_of(path):
    with open(path) as fh:
        for line in fh:
            if line.startswith('DimSize'):
                return tuple(int(x) for x in line.split('=')[1].split())
    return None


def fetch(eid):
    """Download and unpack one grid. Returns None on success, else the reason."""
    last = ''
    for attempt in range(RETRIES + 1):
        try:
            with urllib.request.urlopen(URL.format(eid), timeout=TIMEOUT) as fh:
                blob = fh.read()
            z = zipfile.ZipFile(io.BytesIO(blob))
            names = z.namelist()
            for want, ext in (('energy.mhd', '.mhd'), ('energy.raw', '.raw')):
                src = next((n for n in names if n.endswith(want)), None)
                if src is None:
                    return f'no {want} in the zip'
                with open(os.path.join(DEST, f'{eid}_energy{ext}'), 'wb') as out:
                    out.write(z.read(src))
            dims = dims_of(os.path.join(DEST, f'{eid}_energy.mhd'))
            if dims != GRID_DIMS:
                return f'grid is {dims}, not the reference {GRID_DIMS}'
            return None
        except urllib.error.HTTPError as why:
            return f'HTTP {why.code}'           # a 404 will not improve on retry
        except Exception as why:                # timeouts and truncated reads do
            last = f'{type(why).__name__}: {why}'
            time.sleep(2 * (attempt + 1))
    return last


def main():
    with open(PANEL, newline='', encoding='utf-8') as fh:
        rows = list(csv.DictReader(fh))
    todo = [r for r in rows if not already_there(r['experiment_id'])]
    print(f'{len(rows)} experiments in the panel, {len(rows) - len(todo)} already on disk, '
          f'{len(todo)} to fetch', flush=True)

    failures, done, t0 = [], 0, time.time()
    for i, r in enumerate(todo, 1):
        why = fetch(r['experiment_id'])
        if why:
            failures.append(dict(symbol=r['symbol'], role=r['role'],
                                 experiment_id=r['experiment_id'], plane=r['plane'],
                                 reason=why))
            for ext in ('.mhd', '.raw'):        # never leave half a grid behind
                p = os.path.join(DEST, f'{r["experiment_id"]}_energy{ext}')
                if os.path.exists(p):
                    os.remove(p)
        else:
            done += 1
        if i % 25 == 0 or i == len(todo):
            rate = i / max(time.time() - t0, 1e-9)
            print(f'  {i}/{len(todo)}  {done} ok  {len(failures)} failed  '
                  f'{rate * 60:.0f}/min  eta {(len(todo) - i) / max(rate, 1e-9) / 60:.0f} min',
                  flush=True)

    path = os.path.join(OUT, 'fetch_failures.csv')
    with open(path, 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=['symbol', 'role', 'experiment_id',
                                           'plane', 'reason'])
        w.writeheader(); w.writerows(failures)
    print(f'\n{done} fetched, {len(failures)} failed -> {path}')
    usable = sum(1 for r in rows if already_there(r['experiment_id']))
    print(f'{usable} of {len(rows)} panel experiments are now on disk')


if __name__ == '__main__':
    main()
