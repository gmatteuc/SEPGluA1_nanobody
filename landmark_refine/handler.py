"""
One request in, one response out. Shared by cli.py (one process per call) and
serve.py (one process answering many), so the two can never drift apart.

Request and response are plain dicts of numpy arrays / strings, the shape
scipy.io.loadmat and savemat produce and consume -- see cli.py for the fields.
"""

import os
import sys
import traceback

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np  # noqa: E402


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


def handle(req):
    """Never raises: a failure comes back as ok=0 with the message."""
    resp = dict(ok=0, message='', atlas_pts=np.zeros((0, 2)), hist_pts=np.zeros((0, 2)))
    try:
        import core
        mode = _str(req, 'mode', 'refine').strip().lower()
        atlas, hist = req['atlas'], req['hist']
        labels = req['labels'] if 'labels' in req and np.asarray(req['labels']).size else None

        if mode == 'refine':
            pts = np.asarray(req['atlas_pts'], float).reshape(-1, 2)
            out = core.refine(atlas, hist, pts, labels, snap_r=float(_scalar(req, 'snap_r', core.SNAP_R)))
        elif mode == 'suggest':
            out = core.suggest(atlas, hist, labels, n=int(_scalar(req, 'n', 10)),
                               mirror=bool(_scalar(req, 'mirror', 1)))
        elif mode == 'ping':
            out = dict(ok=True, message='alive')
        else:
            raise ValueError(f'unknown mode {mode!r}')

        for k, v in out.items():
            if isinstance(v, np.ndarray) and v.dtype == bool:
                v = v.astype(float)
            resp[k] = v
        resp['ok'] = int(bool(out['ok']))
        resp['message'] = out.get('message', '')
    except Exception as e:          # noqa: BLE001
        resp['ok'] = 0
        resp['message'] = f'{type(e).__name__}: {e}'
        traceback.print_exc()
    resp['message'] = resp['message'] or ' '
    return resp
