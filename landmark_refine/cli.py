"""
One-shot entry point, for calling from MATLAB via system() when no worker is
running. Pays Python start-up, torch import and model load on every call
(2-3 s); serve.py pays them once.

    python cli.py request.mat response.mat

request.mat (MATLAB -v7):
    mode       'refine' | 'suggest' | 'ping'
    atlas      HxW image (uint8 or float), the atlas plane as shown in the GUI
    hist       HxW image, the histology slice as shown in the GUI (same size)
    atlas_pts  Nx2 [x y], refine only: the points to carry over
    labels     HxW integer label plane (optional; enables boundary snapping
               and boundary scoring)
    snap_r     scalar px (optional, default 12)
    n          scalar, suggest only: how many points (default 10)
    mirror     scalar 0/1, suggest only: mirror pairs (default 1)

response.mat:
    ok, message, atlas_pts, hist_pts, and per-point diagnostics
    (moved, n_local, local_rms for refine; score for suggest),
    n_matches, n_inliers, det, mean_conf, seconds.

Exit code 0 on success, 1 on failure; a response is written either way so
the caller never has to parse stderr.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from scipy.io import loadmat, savemat   # noqa: E402
from handler import handle              # noqa: E402


def main(req_path, resp_path):
    try:
        req = loadmat(req_path)
    except Exception as e:              # noqa: BLE001
        resp = dict(ok=0, message=f'could not read request: {type(e).__name__}: {e}')
    else:
        resp = handle(req)
    savemat(resp_path, resp, do_compression=False)
    return 0 if resp['ok'] else 1


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(2)
    sys.exit(main(sys.argv[1], sys.argv[2]))
