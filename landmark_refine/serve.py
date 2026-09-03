"""
Persistent worker: import torch and load the model ONCE, then answer requests
dropped into a folder. Cuts a GUI proposal from several seconds -- almost all
of it Python start-up -- to the ~0.4 s the matching itself takes.

    python serve.py <folder>

Protocol, all inside <folder>:
    req_<id>.mat    a request, written by the caller under a temporary name
                    and renamed into place, so it is never read half-written
    resp_<id>.mat   the answer, written the same way; the request is deleted
    heartbeat       touched every loop; a caller treats the worker as alive
                    if this is less than a few seconds old
    pid             this process's id, for the record
    stop            create it and the worker exits (and removes it)

The request/response contents are exactly cli.py's. A caller that finds no
live heartbeat falls back to cli.py and gets the same answer, slower.
"""

import os
import sys
import time
import glob

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np                      # noqa: E402
from scipy.io import loadmat, savemat   # noqa: E402
from handler import handle              # noqa: E402

POLL = 0.03          # s between looks at the folder
HEARTBEAT = 1.0      # s between touches


def touch(path):
    with open(path, 'a'):
        os.utime(path, None)


def main(folder):
    os.makedirs(folder, exist_ok=True)
    # Write our own log, unbuffered, rather than trust whoever launched us to
    # have captured stdout: the process may have no console at all (pythonw),
    # and a start-up crash that leaves an empty log is the worst kind.
    log = open(os.path.join(folder, 'worker.log'), 'a', buffering=1)
    sys.stdout = sys.stderr = log
    print(f'--- worker starting, pid {os.getpid()}, python {sys.version.split()[0]}', flush=True)
    try:
        return _serve(folder)
    except BaseException:           # noqa: BLE001  -- the log must say what happened
        import traceback
        traceback.print_exc()
        raise


def _serve(folder):
    stop = os.path.join(folder, 'stop')
    if os.path.exists(stop):
        os.remove(stop)
    with open(os.path.join(folder, 'pid'), 'w') as f:
        f.write(str(os.getpid()))

    # Pay the start-up now, not on the first request: import, model load, and
    # one matching pass so cuDNN has picked its kernels.
    t0 = time.time()
    import core
    core._loftr()
    rng = np.random.default_rng(0)
    im = rng.random((200, 280)).astype(np.float32)
    core.match_multiscale(im, im, sigmas=(0,))
    print(f'landmark_refine worker ready on {core.device()} after {time.time() - t0:.1f} s, '
          f'watching {folder}', flush=True)

    last_beat = 0.0
    while True:
        now = time.time()
        if now - last_beat > HEARTBEAT:
            touch(os.path.join(folder, 'heartbeat'))
            last_beat = now
        if os.path.exists(stop):
            os.remove(stop)
            print('stop requested, exiting', flush=True)
            return 0

        for req_path in sorted(glob.glob(os.path.join(folder, 'req_*.mat'))):
            rid = os.path.basename(req_path)[4:-4]
            t = time.time()
            try:
                req = loadmat(req_path)
                resp = handle(req)
            except Exception as e:          # noqa: BLE001
                resp = dict(ok=0, message=f'{type(e).__name__}: {e}',
                            atlas_pts=np.zeros((0, 2)), hist_pts=np.zeros((0, 2)))
            tmp = os.path.join(folder, f'tmp_resp_{rid}.mat')
            savemat(tmp, resp, do_compression=False)
            os.replace(tmp, os.path.join(folder, f'resp_{rid}.mat'))
            try:
                os.remove(req_path)
            except OSError:
                pass
            print(f'{rid}: {req.get("mode", ["?"])[0] if "mode" in req else "?"} '
                  f'ok={resp["ok"]} in {time.time() - t:.2f} s', flush=True)
        time.sleep(POLL)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(2)
    sys.exit(main(sys.argv[1]))
