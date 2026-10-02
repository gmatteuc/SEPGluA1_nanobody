"""Download the Allen expression grids the ontology panel asks for, once.

One zip per experiment from `api.brain-map.org/grid_data/download/<id>`, holding
`energy.mhd` and `energy.raw`, unpacked into atlas_ish/ under the experiment id,
beside the grids already there. A grid already on disk is left alone, so a re-run
fetches only what is missing, after a partial run too, which matters with several
hundred grids.

Two failures are recorded with their reason in fetch_failures.csv rather than
swallowed: a grid that cannot be downloaded (an HTTP error, not retried: some
experiments have no grid at all), and a grid outside the shared 67 x 41 x 58
reference box, which cannot be placed against the atlas (ish.regions drops such
grids for the same reason). The files of a failed grid are removed, so no half
grid is left, and the panel's real size is a number on disk, not a guess.

Run by run_panel_fetch.py.
"""

import csv
import io
import os
import time
import urllib.error
import urllib.request
import zipfile

from sepmap.config import DATA

PANEL = os.path.join(DATA, "adult_v2", "panel", "panel_v2.csv")
DEST = os.path.join(DATA, "atlas_ish")
OUT = os.path.join(DATA, "adult_v2", "panel")

# one zip per experiment id
URL = "http://api.brain-map.org/grid_data/download/{}"

# the shared reference box, in the header's (x, y, z) order
GRID_DIMS = (67, 41, 58)

# seconds a download may take
TIMEOUT = 180

# further attempts after a timeout or a truncated read
RETRIES = 2


def already_there(eid):
    """Whether both files of the experiment's grid are on disk."""
    return all(
        os.path.exists(os.path.join(DEST, f"{eid}_energy{ext}"))
        for ext in (".mhd", ".raw")
    )


def dims_of(path):
    """The DimSize of a MetaImage header as a tuple, or None when it has none."""
    with open(path) as fh:
        for line in fh:
            if line.startswith("DimSize"):
                return tuple(int(x) for x in line.split("=")[1].split())
    return None


def fetch(eid):
    """Download and unpack one grid. Returns None on success, else the reason."""
    last = ""
    for attempt in range(RETRIES + 1):
        try:
            with urllib.request.urlopen(URL.format(eid), timeout=TIMEOUT) as fh:
                blob = fh.read()
            z = zipfile.ZipFile(io.BytesIO(blob))
            names = z.namelist()
            for want, ext in (("energy.mhd", ".mhd"), ("energy.raw", ".raw")):
                src = next((n for n in names if n.endswith(want)), None)
                if src is None:
                    return f"no {want} in the zip"
                with open(os.path.join(DEST, f"{eid}_energy{ext}"), "wb") as out:
                    out.write(z.read(src))
            dims = dims_of(os.path.join(DEST, f"{eid}_energy.mhd"))
            if dims != GRID_DIMS:
                return f"grid is {dims}, not the reference {GRID_DIMS}"
            return None
        except urllib.error.HTTPError as why:
            # a 404 will not improve on retry
            return f"HTTP {why.code}"
        except Exception as why:
            # timeouts and truncated reads may, so they are retried
            last = f"{type(why).__name__}: {why}"
            time.sleep(2 * (attempt + 1))
    return last


def main():
    """Download the panel's missing grids and write the failures."""
    # the experiments of the panel not yet on disk
    with open(PANEL, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    todo = [r for r in rows if not already_there(r["experiment_id"])]
    print(
        f"{len(rows)} experiments in the panel, {len(rows) - len(todo)} already on disk, "
        f"{len(todo)} to fetch",
        flush=True,
    )

    # one grid after another
    failures, done, t0 = [], 0, time.time()
    for i, r in enumerate(todo, 1):
        why = fetch(r["experiment_id"])
        if why:
            failures.append(
                dict(
                    symbol=r["symbol"],
                    role=r["role"],
                    experiment_id=r["experiment_id"],
                    plane=r["plane"],
                    reason=why,
                )
            )

            # never leave half a grid behind
            for ext in (".mhd", ".raw"):
                p = os.path.join(DEST, f"{r['experiment_id']}_energy{ext}")
                if os.path.exists(p):
                    os.remove(p)
        else:
            done += 1

        # progress every 25 grids
        if i % 25 == 0 or i == len(todo):
            rate = i / max(time.time() - t0, 1e-9)
            print(
                f"  {i}/{len(todo)}  {done} ok  {len(failures)} failed  "
                f"{rate * 60:.0f}/min  "
                f"eta {(len(todo) - i) / max(rate, 1e-9) / 60:.0f} min",
                flush=True,
            )

    # the failures with their reason, and what is on disk now
    path = os.path.join(OUT, "fetch_failures.csv")
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(
            fh, fieldnames=["symbol", "role", "experiment_id", "plane", "reason"]
        )
        w.writeheader()
        w.writerows(failures)
    print(f"\n{done} fetched, {len(failures)} failed -> {path}")
    usable = sum(1 for r in rows if already_there(r["experiment_id"]))
    print(f"{usable} of {len(rows)} panel experiments are now on disk")
