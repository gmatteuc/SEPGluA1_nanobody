"""The automatic annotation, from MATLAB: one process per call.

    python cli.py propose <lightsuite folder>
        the whole brain. Reads, from that folder:
          volume_for_inspection.tiff   the sections, as the GUI shows them
          auto_atlas_planes.mat        the warped atlas the GUI draws (tv, planes x
                                       H x W uint8), written when anchors are saved
          plane_anchors.mat            anchor_slices, anchor_planes (1-based), from
                                       the GUI
        Writes, into the same folder:
          auto_proposal_controlpoints.mat   histology_control_points and
                                            atlas_control_points, one N x 4 cell per
                                            slice as the GUI stores them
                                            ([slice|plane, y, x, t], 1-based)
          auto_proposal_info.mat            planes used, confidence per point,
                                            flags, model version, date
        Neither name ends in 'tform.mat': registerSlicesToAtlas globs for that, and
        a proposal must never be picked up as an annotation.

    python cli.py section <lightsuite folder> <slice> <plane> <response.mat>
        one section at a given plane (both 1-based), for the GUI when the plane of
        a proposed section is changed by hand. Writes atlas, hist (n x 2, 1-based
        y x) and low (n,) to <response.mat>.

    python cli.py sections <lightsuite folder> <request.mat> <response.mat>
        several sections at once (request: slices, planes, 1-based), for the GUI's
        U after an anchor was corrected. Response: atlas, hist, low as cells.

Planes between anchors are interpolated linearly, and extrapolated past the end
ones, as registerSlicesToAtlas does for sections without points. Run by
auto_annotate.m, in the engine's own environment (setup.ps1).
"""

import datetime
import os
import sys
import time

import numpy as np
from scipy.io import loadmat, savemat

# core.py sits beside this script, which MATLAB runs by its path
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import core  # noqa: E402

# the trained models and their VERSION.txt, tracked with the code
WEIGHTS = os.path.join(HERE, "weights")


def read_sections(folder):
    """The sections as the GUI shows them, (S, H, W, 3) uint8."""
    import tifffile

    return tifffile.imread(os.path.join(folder, "volume_for_inspection.tiff"))


def read_atlas(folder):
    """The warped atlas the GUI saved (auto_atlas_planes.mat), planes x H x W uint8."""
    f = os.path.join(folder, "auto_atlas_planes.mat")
    if not os.path.exists(f):
        raise SystemExit(f"no {f}: open the brain in the GUI and save the anchors first")
    try:
        return loadmat(f, variable_names=["tv"])["tv"]
    except NotImplementedError:
        # a -v7.3 file, which loadmat cannot read: HDF5, with the axes reversed
        import h5py

        with h5py.File(f, "r") as fh:
            return np.array(fh["tv"]).transpose()


def planes_from_anchors(folder, n_slices, n_planes):
    """The atlas plane of every slice, and the anchor slices and planes, all 1-based.

    The planes are read from plane_anchors.mat when the GUI saved them for every
    slice; otherwise they are interpolated between the anchors and extended past
    the end ones. Each is clipped to the atlas's 1 to `n_planes`.
    """
    a = loadmat(os.path.join(folder, "plane_anchors.mat"), squeeze_me=True)
    s = np.atleast_1d(a["anchor_slices"]).astype(float)
    p = np.atleast_1d(a["anchor_planes"]).astype(float)

    # the GUI saves the plane its own interpolation gives every slice: use those,
    # so what is proposed is what was on screen
    if "planes" in a and np.size(a["planes"]) == n_slices:
        planes = np.clip(
            np.rint(np.atleast_1d(a["planes"]).astype(float)), 1, n_planes
        ).astype(int)
        return planes, s.astype(int), p.astype(int)
    if len(s) < 2:
        raise SystemExit("at least two anchor sections are needed")
    order = np.argsort(s)
    s, p = s[order], p[order]
    x = np.arange(1, n_slices + 1, dtype=float)

    # linear inside, extended past the ends along the nearest pair of anchors
    y = np.interp(x, s, p)
    lo, hi = x < s[0], x > s[-1]
    y[lo] = p[0] + (x[lo] - s[0]) * (p[1] - p[0]) / (s[1] - s[0])
    y[hi] = p[-1] + (x[hi] - s[-1]) * (p[-1] - p[-2]) / (s[-1] - s[-2])
    return np.clip(np.rint(y), 1, n_planes).astype(int), s.astype(int), p.astype(int)


def model_version():
    """The text of weights/VERSION.txt, or "unversioned" without one."""
    f = os.path.join(WEIGHTS, "VERSION.txt")
    return (
        open(f, encoding="utf-8").read().strip() if os.path.exists(f) else "unversioned"
    )


def propose(folder):
    """Propose every slice of the brain in `folder` and write the two proposal files."""
    t0 = time.time()

    # the sections, the atlas planes they sit on, and the proposal
    vol = read_sections(folder)
    tv = read_atlas(folder)
    n = vol.shape[0]
    planes, anchor_s, anchor_p = planes_from_anchors(folder, n, tv.shape[0])
    lnet, mnet = core.load_models(
        os.path.join(WEIGHTS, "landmark.pt"), os.path.join(WEIGHTS, "matcher.pt")
    )
    res = core.propose(
        [vol[k] for k in range(n)], [tv[planes[k] - 1] for k in range(n)], lnet, mnet
    )

    # one N x 4 cell per slice, as the GUI stores points: [slice or plane, y, x, t],
    # 1-based
    hist = np.empty((n, 1), object)
    atl = np.empty((n, 1), object)
    spread = np.empty((n, 1), object)
    low = np.empty((n, 1), object)
    for k, r in enumerate(res):
        m = len(r["atlas"])
        hist[k, 0] = np.column_stack([np.full(m, k + 1), r["hist"] + 1, np.zeros(m)])
        atl[k, 0] = np.column_stack([np.full(m, planes[k]), r["atlas"] + 1, np.zeros(m)])
        spread[k, 0] = r["spread"].reshape(-1, 1)
        low[k, 0] = r["low"].reshape(-1, 1).astype(np.uint8)
    savemat(
        os.path.join(folder, "auto_proposal_controlpoints.mat"),
        {"histology_control_points": hist, "atlas_control_points": atl},
        do_compression=True,
    )
    savemat(
        os.path.join(folder, "auto_proposal_info.mat"),
        {
            "planes": planes.reshape(-1, 1),
            "anchor_slices": anchor_s.reshape(-1, 1),
            "anchor_planes": anchor_p.reshape(-1, 1),
            "spread": spread,
            "low_confidence": low,
            "model_version": model_version(),
            "device": core.device(),
            "created": datetime.datetime.now().isoformat(timespec="seconds"),
        },
        do_compression=True,
    )
    print(
        f"auto_annotation: {n} sections, {sum(len(r['atlas']) for r in res)} pairs, "
        f"{core.device()}, {time.time() - t0:.0f} s -> {folder}"
    )


def section(folder, slice_1b, plane_1b, response):
    """Propose one slice at one plane (both 1-based), the answer into `response`."""
    vol = read_sections(folder)
    tv = read_atlas(folder)
    lnet, mnet = core.load_models(
        os.path.join(WEIGHTS, "landmark.pt"), os.path.join(WEIGHTS, "matcher.pt")
    )
    r = core.propose([vol[slice_1b - 1]], [tv[plane_1b - 1]], lnet, mnet)[0]
    savemat(
        response,
        {
            "ok": 1,
            "message": "",
            "atlas": r["atlas"] + 1,
            "hist": r["hist"] + 1,
            "low": r["low"].astype(np.uint8),
        },
    )


def sections(folder, request, response):
    """Propose several slices, each at its plane (the GUI's U), into `response`.

    `request` holds the slices and the planes, 1-based; the answer holds one cell
    per slice.
    """
    q = loadmat(request, squeeze_me=True)
    slices = np.atleast_1d(q["slices"]).astype(int)
    planes = np.atleast_1d(q["planes"]).astype(int)
    vol = read_sections(folder)
    tv = read_atlas(folder)
    lnet, mnet = core.load_models(
        os.path.join(WEIGHTS, "landmark.pt"), os.path.join(WEIGHTS, "matcher.pt")
    )
    res = core.propose(
        [vol[s - 1] for s in slices], [tv[p - 1] for p in planes], lnet, mnet
    )
    atl = np.empty((len(res), 1), object)
    hist = np.empty((len(res), 1), object)
    low = np.empty((len(res), 1), object)
    for k, r in enumerate(res):
        atl[k, 0] = r["atlas"] + 1
        hist[k, 0] = r["hist"] + 1
        low[k, 0] = r["low"].astype(np.uint8).reshape(-1, 1)
    savemat(response, {"ok": 1, "message": "", "atlas": atl, "hist": hist, "low": low})


def main(argv):
    """Run the mode `argv` names; any other call exits with the usage.

    In the section modes an exception goes into the response file, ok 0 and its
    message, which the GUI shows; a SystemExit (no atlas file, too few anchors)
    still ends the process without one.
    """
    if len(argv) >= 2 and argv[0] == "propose":
        propose(argv[1])
    elif len(argv) >= 4 and argv[0] == "sections":
        try:
            sections(argv[1], argv[2], argv[3])
        except Exception as err:
            # the GUI shows the message
            savemat(argv[3], {"ok": 0, "message": f"{type(err).__name__}: {err}"})
    elif len(argv) >= 5 and argv[0] == "section":
        try:
            section(argv[1], int(argv[2]), int(argv[3]), argv[4])
        except Exception as err:
            # the GUI shows the message
            savemat(argv[4], {"ok": 0, "message": f"{type(err).__name__}: {err}"})
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main(sys.argv[1:])
