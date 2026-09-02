"""
Command-line entry point, for calling from MATLAB via system().

    python cli.py request.mat response.mat

request.mat (MATLAB -v7):
    mode       'refine' | 'suggest'
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
    n_matches, n_inliers, det, seconds.

Exit code 0 on success, 1 on failure; a response is written either way so
the caller never has to parse stderr.
"""

import os
import sys
import traceback

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np                      # noqa: E402
from scipy.io import loadmat, savemat   # noqa: E402


def _scalar(d, key, default):
    if key not in d:
        return default
    v = np.asarray(d[key]).ravel()
    return default if v.size == 0 else v[0]


def _str(d, key, default):
    if key not in d:
        return default
    v = d[key]
    return str(v[0]) if isinstance(v, np.ndarray) and v.size else str(v)


def main(req_path, resp_path):
    resp = dict(ok=0, message='', atlas_pts=np.zeros((0, 2)), hist_pts=np.zeros((0, 2)))
    try:
        req = loadmat(req_path)
        mode = _str(req, 'mode', 'refine').strip().lower()
        atlas, hist = req['atlas'], req['hist']
        labels = req['labels'] if 'labels' in req and np.asarray(req['labels']).size else None

        import core
        if mode == 'refine':
            pts = np.asarray(req['atlas_pts'], float).reshape(-1, 2)
            out = core.refine(atlas, hist, pts, labels, snap_r=float(_scalar(req, 'snap_r', core.SNAP_R)))
        elif mode == 'suggest':
            out = core.suggest(atlas, hist, labels, n=int(_scalar(req, 'n', 10)),
                               mirror=bool(_scalar(req, 'mirror', 1)))
        else:
            raise ValueError(f'unknown mode {mode!r}')

        for k, v in out.items():
            if isinstance(v, np.ndarray):
                v = v.astype(float) if v.dtype == bool else v
            resp[k] = v
        resp['ok'] = int(bool(out['ok']))
        resp['message'] = out.get('message', '')
    except Exception as e:          # noqa: BLE001
        resp['ok'] = 0
        resp['message'] = f'{type(e).__name__}: {e}'
        traceback.print_exc()
    resp['message'] = resp['message'] or ' '
    savemat(resp_path, resp, do_compression=False)
    return 0 if resp['ok'] else 1


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(2)
    sys.exit(main(sys.argv[1], sys.argv[2]))
