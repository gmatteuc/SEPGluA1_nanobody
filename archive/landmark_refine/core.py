"""
Landmark refinement for the control-point GUI.

Two operations, both built on the same machinery and both validated against
19 hand-annotated slices of MG903 (P20) vs DeMBA P20 -- see
D:\\sep_histology\\sandbox_landmark_matching\\README.md for every number.

  refine(atlas, hist, atlas_pts, labels)
      The "adjusted carry-over". Keeps the annotator's choice of landmarks
      (the atlas points of the previous slice) and fixes their coordinates:
      the atlas side is nudged onto the nearest structure boundary of the new
      plane, the histology side is discarded and re-derived by pushing each
      atlas point through a local affine fitted to the nearest image matches.
      Measured: 9.7 px median vs 15.2 for a plain copy, beats it on 16/18
      slices, keeps the annotator's mirror-pair structure.

  suggest(atlas, hist, labels, n, mirror)
      Points from scratch, for a slice with nothing to carry. Multi-scale
      matches, scored by confidence + boundary proximity + cornerness, chosen
      as left/right pairs across medial/mid/lateral bands.
      Measured: 10.2 px median at 10 points, 9.3 at 15 unpaired.

Matching is LoFTR (kornia, pretrained 'outdoor', nothing fine-tuned), run at
four blur levels and pooled -- pooling was the single largest improvement
found. The affine RANSAC refuses reflections: a mirrored solution fits a
near-symmetric brain perfectly and is a valid affine, so it has to be excluded
by the model, not by the residual.

Coordinates everywhere here are (x, y) in pixels of the images as given. The
MATLAB wrapper handles the GUI's [plane y x t] convention.
"""

import os
import time

import numpy as np
import torch
from scipy.ndimage import distance_transform_edt, gaussian_filter, sobel, uniform_filter

_LOFTR = None
_DEVICE = None
SIGMAS = (0, 1, 2, 3)          # suggest: measured at 9.3 px with all four
SIGMAS_REFINE = (0, 2)         # refine: (0,2) measured 9.5 px vs 9.7 for all four, at half the cost
W_BOUND = 3.0          # boundary bonus; chosen by leave-one-slice-out CV
W_CORNER = 1.0
SNAP_R = 12            # px; snapping is a nudge, accuracy comes from the transfer
K_LOCAL = 12           # matches per local affine
LOCAL_RADIUS = 60.0    # px
MIRROR_TOL = 14.0      # px, partner must sit this close to the mirrored position
BANDS = (0, 45, 110, 10 ** 6)


# ------------------------------------------------------------------ images

def as_float(im):
    im = np.asarray(im)
    if im.ndim == 3:
        im = im[..., 0]
    im = im.astype(np.float32)
    if im.max() > 1.0:
        im = im / 255.0
    return np.clip(im, 0, 1)


def _pad8(im):
    h, w = im.shape
    H, W = (h + 7) // 8 * 8, (w + 7) // 8 * 8
    out = np.zeros((H, W), np.float32)
    out[:h, :w] = im
    return out


def device():
    """The GPU when there is one, else the CPU. Four LoFTR passes take ~2 s on
    CPU and a fraction of that on a modest GPU; nothing else here is slow."""
    global _DEVICE
    if _DEVICE is None:
        _DEVICE = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    return _DEVICE


WEIGHTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'weights', 'loftr_outdoor.ckpt')


def _loftr():
    """The model, loaded once. Weights come from the bundled file next to this
    module, so a machine with no internet works; only if that file is missing
    does kornia go and download them."""
    global _LOFTR
    if _LOFTR is None:
        import kornia.feature as KF
        torch.set_grad_enabled(False)
        if os.path.exists(WEIGHTS):
            m = KF.LoFTR(pretrained=None)
            m.load_state_dict(torch.load(WEIGHTS, map_location='cpu')['state_dict'])
        else:
            m = KF.LoFTR(pretrained='outdoor')
        _LOFTR = m.eval().to(device())
    return _LOFTR


def match_multiscale(atlas, hist, sigmas=SIGMAS):
    """Pooled LoFTR matches across blur levels. Returns (pa, ph, conf)."""
    m = _loftr()
    dev = device()
    pas, phs, cfs = [], [], []
    for s in sigmas:
        a = gaussian_filter(atlas, s) if s else atlas
        b = gaussian_filter(hist, s) if s else hist
        ta = torch.from_numpy(_pad8(a))[None, None].to(dev)
        tb = torch.from_numpy(_pad8(b))[None, None].to(dev)
        out = m({'image0': ta, 'image1': tb})
        pas.append(out['keypoints0'].cpu().numpy())
        phs.append(out['keypoints1'].cpu().numpy())
        cfs.append(out['confidence'].cpu().numpy())
    return np.vstack(pas), np.vstack(phs), np.concatenate(cfs)


# ------------------------------------------------------------------ affine

def fit_affine(src, dst, w=None):
    n = len(src)
    A = np.zeros((2 * n, 6))
    A[0::2, 0:2] = src; A[0::2, 2] = 1
    A[1::2, 3:5] = src; A[1::2, 5] = 1
    b = np.asarray(dst, float).reshape(-1)
    if w is not None:
        ww = np.repeat(np.sqrt(w), 2)
        A, b = A * ww[:, None], b * ww
    x, *_ = np.linalg.lstsq(A, b, rcond=None)
    return np.array([[x[0], x[1], x[2]], [x[3], x[4], x[5]]])


def apply_affine(M, pts):
    return np.asarray(pts, float) @ M[:, :2].T + M[:, 2]


def ransac_affine(src, dst, thresh=6.0, iters=3000, seed=0):
    """Affine RANSAC that refuses reflections (det <= 0)."""
    rng = np.random.default_rng(seed)
    n = len(src)
    if n < 3:
        return None, np.zeros(n, bool)
    best_in, best_M = np.zeros(n, bool), None
    for _ in range(iters):
        idx = rng.choice(n, 3, replace=False)
        try:
            M = fit_affine(src[idx], dst[idx])
        except np.linalg.LinAlgError:
            continue
        det = np.linalg.det(M[:, :2])
        if det <= 0 or det < 0.05 or det > 20:
            continue
        inl = np.linalg.norm(apply_affine(M, src) - dst, axis=1) < thresh
        if inl.sum() > best_in.sum():
            best_in, best_M = inl, M
    if best_M is None or best_in.sum() < 3:
        return None, np.zeros(n, bool)
    return fit_affine(src[best_in], dst[best_in]), best_in


# ----------------------------------------------------------------- anatomy

def boundary_distance(labels):
    """Distance to the nearest internal structure boundary, from a label plane."""
    lab = np.asarray(labels).astype(np.int64)
    inside = lab > 0
    diff = np.zeros_like(inside)
    for sh in ((0, 1), (0, -1), (1, 0), (-1, 0)):
        r = np.roll(lab, sh, axis=(0, 1))
        diff |= (r != lab) & inside & (r > 0)
    return distance_transform_edt(~diff)


def cornerness(im, win=7):
    """Shi-Tomasi: smaller structure-tensor eigenvalue. High at corners, low on edges."""
    gx = sobel(im, axis=1, mode='nearest')
    gy = sobel(im, axis=0, mode='nearest')
    Axx, Ayy, Axy = (uniform_filter(gx * gx, win), uniform_filter(gy * gy, win),
                     uniform_filter(gx * gy, win))
    tr, det = Axx + Ayy, Axx * Ayy - Axy * Axy
    return tr / 2 - np.sqrt(np.maximum(0.0, tr * tr / 4 - det))


def find_midline(im):
    """Column that best maps the image onto its own mirror image."""
    h, w = im.shape
    col = im - im.mean()
    best, best_x = -np.inf, w / 2
    for x in range(int(0.35 * w), int(0.65 * w)):
        half = min(x, w - x)
        if half < 40:
            continue
        L, R = col[:, x - half:x], col[:, x:x + half][:, ::-1]
        d = np.linalg.norm(L) * np.linalg.norm(R)
        if d == 0:
            continue
        s = float((L * R).sum() / d)
        if s > best:
            best, best_x = s, x
    return best_x


def _sample(img, pts):
    y = np.clip(np.round(pts[:, 1]).astype(int), 0, img.shape[0] - 1)
    x = np.clip(np.round(pts[:, 0]).astype(int), 0, img.shape[1] - 1)
    return img[y, x]


# --------------------------------------------------------------- refinement

def snap_to_boundary(pts, dmap, r=SNAP_R):
    """Pull each point to the nearest boundary pixel within r. Returns (pts, moved)."""
    ys, xs = np.nonzero(dmap <= 0.5)
    B = np.column_stack([xs, ys]).astype(float)
    out = np.asarray(pts, float).copy()
    moved = np.zeros(len(pts), bool)
    if len(B) == 0:
        return out, moved
    for i, p in enumerate(out):
        d = np.linalg.norm(B - p, axis=1)
        j = int(np.argmin(d))
        if d[j] <= r:
            out[i], moved[i] = B[j], True
    return out, moved


def transfer(pts_atlas, pa, ph, M_global, k=K_LOCAL, radius=LOCAL_RADIUS):
    """Atlas -> histology through a local affine on the nearest inlier matches.

    Returns (hist_pts, n_local, local_rms): how many matches supported each
    point and how well they agreed, so the GUI can tell a confident transfer
    from a fallback to the global affine.
    """
    out = np.zeros((len(pts_atlas), 2))
    n_local = np.zeros(len(pts_atlas), int)
    rms = np.zeros(len(pts_atlas))
    for i, p in enumerate(np.asarray(pts_atlas, float)):
        d = np.linalg.norm(pa - p, axis=1)
        idx = np.argsort(d)[:k]
        idx = idx[d[idx] <= radius]
        if len(idx) >= 4:
            w = np.exp(-(d[idx] / (radius / 2)) ** 2)
            M = fit_affine(pa[idx], ph[idx], w)
            rms[i] = float(np.sqrt(np.mean(np.sum((apply_affine(M, pa[idx]) - ph[idx]) ** 2, 1))))
        else:
            M = M_global
            rms[i] = np.nan
        n_local[i] = len(idx)
        out[i] = apply_affine(M, p[None])[0]
    return out, n_local, rms


# ---------------------------------------------------------------- selection

def spread_select(pts, score, k, imshape):
    """k high-scoring points that also cover the image (greedy farthest-point)."""
    if len(pts) <= k:
        return np.arange(len(pts))
    order = np.argsort(-score)
    cand = list(order[:max(k * 12, 150)])
    chosen = [cand[0]]
    d_min = 0.12 * max(imshape)
    while len(chosen) < k and d_min > 1.0:
        for i in cand:
            if len(chosen) >= k:
                break
            if i in chosen:
                continue
            if np.linalg.norm(pts[chosen] - pts[i], axis=1).min() > d_min:
                chosen.append(i)
        d_min *= 0.6
    return np.array(chosen[:k])


def select_pairs_stratified(pa, score, xmid, k, tol=MIRROR_TOL):
    """Mirror pairs drawn round-robin across medial / mid / lateral bands.

    Best-first pairing alone drifts to the midline (highest contrast, partner
    always available). A midline point also mirrors onto itself, so partners
    must sit on the opposite side.
    """
    mir = np.column_stack([2 * xmid - pa[:, 0], pa[:, 1]])
    side = np.sign(pa[:, 0] - xmid)
    lat = np.abs(pa[:, 0] - xmid)
    pairs = []
    for i in np.argsort(-score):
        if lat[i] < 5 or side[i] <= 0:
            continue
        d = np.linalg.norm(pa - mir[i], axis=1)
        d[side == side[i]] = np.inf
        d[i] = np.inf
        j = int(np.argmin(d))
        if d[j] <= tol:
            pairs.append((score[i] + score[j], i, j, lat[i]))
    pairs.sort(reverse=True)
    band = lambda l: int(np.searchsorted(BANDS, l, side='right')) - 1
    by_band = {b: [p for p in pairs if band(p[3]) == b] for b in range(len(BANDS) - 1)}
    chosen, anchors, used = [], [], set()
    d_min = 28.0
    while len(chosen) < k:
        progressed = False
        for b in (2, 1, 0):
            if len(chosen) >= k:
                break
            for sc, i, j, l in by_band[b]:
                if i in used or j in used:
                    continue
                if anchors and min(np.linalg.norm(pa[i] - a) for a in anchors) < d_min:
                    continue
                chosen += [i, j]; anchors += [pa[i], pa[j]]; used |= {i, j}
                progressed = True
                break
        if not progressed:
            d_min *= 0.6
            if d_min < 3:
                break
    return np.array(chosen[:k], int)


# ------------------------------------------------------------- public API

MIN_INLIERS = 30       # pooled inliers on real slices were 400-700; blank or wrong images give a handful


def _matches(atlas, hist, sigmas=SIGMAS):
    info = dict(n_matches=0, n_inliers=0, det=0.0, mean_conf=0.0)
    # A blank or constant panel matches itself perfectly everywhere -- LoFTR
    # returns hundreds of identical-feature matches consistent with the
    # identity, which no inlier count would reject. Check the input first.
    if atlas.std() < 1e-3 or hist.std() < 1e-3:
        return None, None, None, None, info
    pa, ph, cf = match_multiscale(atlas, hist, sigmas)
    M, inl = ransac_affine(pa, ph)
    info.update(n_matches=int(len(pa)), n_inliers=int(inl.sum()),
                det=float(np.linalg.det(M[:, :2])) if M is not None else 0.0,
                mean_conf=float(cf[inl].mean()) if inl.any() else 0.0)
    # An affine through three matches always exists; that is not the same as
    # the two images having anything to do with each other. Refuse rather than
    # hand back confident nonsense.
    if M is None or inl.sum() < MIN_INLIERS:
        return None, None, None, None, info
    return pa[inl], ph[inl], cf[inl], M, info


DISAGREE_PX = 15.0     # the two transfers further apart than this: flag the point


def refine(atlas, hist, atlas_pts, labels=None, snap_r=SNAP_R, hist_prev=None, hist_pts_prev=None):
    """Adjusted carry-over. See module docstring.

    With hist_prev and hist_pts_prev -- the previous slice's histology image
    and the annotator's histology points on it -- a second, independent
    transfer is made by matching histology to HISTOLOGY (same modality,
    adjacent sections: far easier than template to tissue) and the two are
    averaged. Measured per landmark against the annotator's final positions:
    plain copy 12.5 px of correction, atlas->hist alone 9.7, the average 7.2,
    better on 28 of 33 slices. Where the two transfers disagree the annotator
    had to move the point more (r = +0.51), so `disagree` is returned per point
    and `uncertain` marks those beyond DISAGREE_PX.

    Returns dict: atlas_pts (Nx2), hist_pts (Nx2), moved (N bool), n_local (N),
    local_rms (N), disagree (N), uncertain (N bool), ok, message, match info.
    """
    t0 = time.time()
    atlas, hist = as_float(atlas), as_float(hist)
    atlas_pts = np.asarray(atlas_pts, float).reshape(-1, 2)
    n = len(atlas_pts)
    empty = dict(atlas_pts=atlas_pts, hist_pts=np.full_like(atlas_pts, np.nan),
                 moved=np.zeros(n, bool), n_local=np.zeros(n, int), local_rms=np.full(n, np.nan),
                 disagree=np.full(n, np.nan), uncertain=np.zeros(n, bool))
    pa, ph, cf, M, info = _matches(atlas, hist, SIGMAS_REFINE)
    if M is None:
        return dict(ok=False, message='no affine found: too few consistent matches',
                    seconds=time.time() - t0, **empty, **info)
    if labels is not None and snap_r > 0:
        atlas_pts, moved = snap_to_boundary(atlas_pts, boundary_distance(labels), snap_r)
    else:
        moved = np.zeros(n, bool)
    hist_pts, n_local, rms = transfer(atlas_pts, pa, ph, M)

    disagree = np.full(n, np.nan)
    info['n_inliers_hh'] = 0
    if hist_prev is not None and hist_pts_prev is not None and len(hist_pts_prev) == n:
        qa, qh, _, Mh, info_hh = _matches(as_float(hist_prev), hist, SIGMAS_REFINE)
        info['n_inliers_hh'] = info_hh['n_inliers']
        if Mh is not None:
            h2, _, _ = transfer(np.asarray(hist_pts_prev, float).reshape(-1, 2), qa, qh, Mh)
            disagree = np.linalg.norm(hist_pts - h2, axis=1)
            hist_pts = 0.5 * (hist_pts + h2)
    uncertain = np.nan_to_num(disagree, nan=0.0) > DISAGREE_PX

    return dict(ok=True, message='', atlas_pts=atlas_pts, hist_pts=hist_pts, moved=moved,
                n_local=n_local, local_rms=rms, disagree=disagree, uncertain=uncertain,
                seconds=time.time() - t0, **info)


def suggest(atlas, hist, labels=None, n=10, mirror=True):
    """Points from scratch. See module docstring.

    Returns dict: atlas_pts (Mx2), hist_pts (Mx2), score (M), ok, message, info.
    """
    t0 = time.time()
    atlas, hist = as_float(atlas), as_float(hist)
    pa, ph, cf, M, info = _matches(atlas, hist)
    if M is None:
        return dict(ok=False, message='no affine found: too few consistent matches',
                    atlas_pts=np.zeros((0, 2)), hist_pts=np.zeros((0, 2)), score=np.zeros(0),
                    seconds=time.time() - t0, **info)
    score = (cf - cf.mean()) / (cf.std() + 1e-9)
    if labels is not None:
        score = score + W_BOUND * np.exp(-_sample(boundary_distance(labels), pa) / 4.0)
    ca = _sample(cornerness(atlas), pa)
    score = score + W_CORNER * ca / (ca.max() + 1e-12)
    if mirror:
        sel = select_pairs_stratified(pa, score, find_midline(atlas), n)
        if len(sel) < min(n, 6):                    # pairing starved: fall back
            sel = spread_select(pa, score, n, atlas.shape)
    else:
        sel = spread_select(pa, score, n, atlas.shape)
    return dict(ok=True, message='', atlas_pts=pa[sel], hist_pts=ph[sel], score=score[sel],
                seconds=time.time() - t0, **info)
