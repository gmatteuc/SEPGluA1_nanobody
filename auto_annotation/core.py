"""
The automatic control-point annotation: the algorithms, no file I/O.

Ported from the validation sandbox (SEPGluA1_autoannotation) with every setting
unchanged; the numbers quoted here are that sandbox's, measured on the 24
hand-annotated brains with models that never saw the brain being scored.

For each section, given the atlas plane it sits on:

  1  map       image-only registration of the section onto its atlas plane:
               an affine, then a smooth deformation on a 12 x 17 control grid,
               both by normalised gradient fields (edge alignment) on the DAPI
               channel. No points involved. Lands each point at ~8 px, the level
               of an affine fitted to a human's own other points.
  2  where     the landmark network (a small U-Net) reads the atlas plane alone
               and paints where a human would click; its 32 strongest peaks, at
               least 20 px apart, are the landmarks.
  3  correct   the matcher (a small ConvNet) compares a 96 px atlas patch with
               the registered-section patch at each landmark and predicts the
               shift the registration missed; averaged over 9 shifted crops.
               8.1 -> 7.15 px on held-out brains.
  4  pair      partner point on the section = registration(landmark + shift).

Each pair also gets a confidence: how much the 9 shifted answers disagree
(px). The least confident quarter of each section is flagged for review.

Coordinates are (y, x), 0-based, in the 400 x 570 frame the GUI shows.
"""

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

# ------------------------------------------------------------------ settings

# registration (sandbox register_affine.py / register_deform.py / full_auto.py)
AFFINE_SCALES = ((4, 150, 3e-2), (2, 100, 1e-2), (1, 60, 5e-3))   # (downsample, steps, learning rate)
BLUR = {4: 1.0, 2: 1.5, 1: 2.0}
GRID = (12, 17)                       # deformation control points, rows x cols
DEFORM_SCALES = ((2, 120, 2e-3), (1, 80, 1e-3))
BENDING = 100.0                       # stiffness of the deformation
BATCH = 48                            # sections registered together on the GPU

# landmarks (sandbox landmark_net.py)
N_POINTS = 32
MIN_DIST = 20                         # px between landmarks
TISSUE = 12                           # atlas grey level (0-255) above which is brain

# matcher (sandbox learn_match.py)
P = 96                                # patch the network sees
SHIFT = 12                            # stored margin around it
STORE = P + 2 * SHIFT
TTA = [(dy, dx) for dy in (-6, 0, 6) for dx in (-6, 0, 6)]

LOW_CONFIDENCE = 0.25                 # least confident share of each section, flagged


def device():
    return 'cuda' if torch.cuda.is_available() else 'cpu'


# -------------------------------------------------------------------- images

def section_image(rgb):
    """A volume_for_inspection page (H, W, 3) -> the DAPI channel, 0-1 (1st-99.5th percentile)."""
    f = rgb.astype(np.float32)[..., 2]
    lo, hi = np.percentile(f, 1), np.percentile(f, 99.5)
    return np.clip((f - lo) / max(hi - lo, 1e-6), 0, 1).astype(np.float32)


def blur(x, sigma):
    if sigma <= 0:
        return x
    r = int(3 * sigma)
    k = torch.exp(-torch.arange(-r, r + 1, device=x.device, dtype=x.dtype) ** 2 / (2 * sigma ** 2))
    k = k / k.sum()
    x = F.conv2d(x, k.view(1, 1, 1, -1), padding=(0, r))
    return F.conv2d(x, k.view(1, 1, -1, 1), padding=(r, 0))


# -------------------------------------------------------------- registration

def ngf(a, b, m):
    """Normalised gradient fields: do the edges line up, whatever their contrast."""
    def grads(x):
        return x[..., 1:, :-1] - x[..., :-1, :-1], x[..., :-1, 1:] - x[..., :-1, :-1]
    ay, ax = grads(a); by, bx = grads(b)
    mm = m[..., :-1, :-1]

    def eps(gy, gx):
        with torch.no_grad():
            mag = torch.sqrt(gy * gy + gx * gx + 1e-12)
            e = 0.5 * (mag * mm).sum(dim=(1, 2, 3), keepdim=True) / mm.sum(dim=(1, 2, 3), keepdim=True).clamp_min(1)
        return e.clamp_min(1e-4)
    ea, eb = eps(ay, ax), eps(by, bx)
    na = torch.sqrt(ay * ay + ax * ax + ea * ea)
    nb = torch.sqrt(by * by + bx * bx + eb * eb)
    dot = (ay * by + ax * bx) / (na * nb)
    return (dot * dot * mm).sum(dim=(1, 2, 3)) / mm.sum(dim=(1, 2, 3)).clamp_min(1)


def register_affine(moving, fixed):
    """Affine theta (B, 2, 3), normalised coordinates, atlas -> section."""
    dev = moving.device
    b = moving.shape[0]
    eye = torch.tensor([[1., 0., 0.], [0., 1., 0.]], device=dev).expand(b, 2, 3)
    delta = torch.zeros(b, 2, 3, device=dev, requires_grad=True)
    for factor, steps, lr in AFFINE_SCALES:
        mv = F.avg_pool2d(moving, factor) if factor > 1 else moving
        fx = F.avg_pool2d(fixed, factor) if factor > 1 else fixed
        mv, fx = blur(mv, BLUR[factor]), blur(fx, BLUR[factor])
        mask = (fx > 0.05).float()
        opt = torch.optim.Adam([delta], lr=lr)
        for _ in range(steps):
            grid = F.affine_grid(eye + delta, fx.shape, align_corners=False)
            warped = F.grid_sample(mv, grid, align_corners=False, padding_mode='zeros')
            loss = (1 - ngf(warped, fx, mask)).sum()
            opt.zero_grad(); loss.backward(); opt.step()
    return (eye + delta).detach()


def bending(field):
    dyy = field[..., 2:, :] - 2 * field[..., 1:-1, :] + field[..., :-2, :]
    dxx = field[..., :, 2:] - 2 * field[..., :, 1:-1] + field[..., :, :-2]
    return (dyy ** 2).mean(dim=(1, 2, 3)) + (dxx ** 2).mean(dim=(1, 2, 3))


def upsample(field, h, w):
    return F.interpolate(field, size=(h, w), mode='bicubic', align_corners=True)


def register_deform(moving, fixed, theta):
    """A smooth displacement (B, 2, GRID), normalised units, added to the affine."""
    b = moving.shape[0]
    field = torch.zeros(b, 2, *GRID, device=moving.device, requires_grad=True)
    for factor, steps, lr in DEFORM_SCALES:
        mv = F.avg_pool2d(moving, factor) if factor > 1 else moving
        fx = F.avg_pool2d(fixed, factor) if factor > 1 else fixed
        mv, fx = blur(mv, BLUR[factor]), blur(fx, BLUR[factor])
        mask = (fx > 0.05).float()
        base = F.affine_grid(theta, fx.shape, align_corners=False)
        opt = torch.optim.Adam([field], lr=lr)
        for _ in range(steps):
            disp = upsample(field, fx.shape[2], fx.shape[3]).permute(0, 2, 3, 1)
            warped = F.grid_sample(mv, base + disp, align_corners=False, padding_mode='zeros')
            loss = ((1 - ngf(warped, fx, mask)) + BENDING * bending(field)).sum()
            opt.zero_grad(); loss.backward(); opt.step()
    return field.detach()


def register(sections, planes):
    """Every section onto its atlas plane. sections, planes: lists of (H, W) float 0-1."""
    dev = device()
    thetas, fields = [], []
    for s in range(0, len(sections), BATCH):
        mv = torch.tensor(np.stack(sections[s:s + BATCH]), device=dev)[:, None]
        fx = torch.tensor(np.stack(planes[s:s + BATCH]), device=dev)[:, None]
        th = register_affine(mv, fx)
        thetas.append(th)
        fields.append(register_deform(mv, fx, th))
    return torch.cat(thetas), torch.cat(fields)


def warp_to_atlas(section, theta, field):
    """The section resampled into the atlas frame by the registration."""
    h, w = section.shape
    grid = F.affine_grid(theta[None], (1, 1, h, w), align_corners=False) + \
        upsample(field[None], h, w).permute(0, 2, 3, 1)
    img = torch.tensor(section, device=theta.device)[None, None]
    return F.grid_sample(img, grid, align_corners=False)[0, 0].cpu().numpy()


def map_points(theta, field, pts, h, w):
    """Atlas (y, x) points through affine + displacement into the section."""
    theta = theta.cpu().numpy() if torch.is_tensor(theta) else theta
    xn = (2 * pts[:, 1] + 1) / w - 1
    yn = (2 * pts[:, 0] + 1) / h - 1
    xa, ya = theta @ np.vstack([xn, yn, np.ones_like(xn)])
    full = upsample(field[None], h, w)
    where = torch.tensor(np.column_stack([xn, yn]), dtype=torch.float32, device=field.device).view(1, 1, -1, 2)
    d = F.grid_sample(full, where, align_corners=False)[0, :, 0].cpu().numpy()
    return np.column_stack([((ya + d[1] + 1) * h - 1) / 2, ((xa + d[0] + 1) * w - 1) / 2])


# ----------------------------------------------------------------- landmarks

class LandmarkNet(nn.Module):
    """Small U-Net: atlas plane -> where a human would click (logits)."""

    def __init__(self, c=(16, 32, 64, 96)):
        super().__init__()

        def block(i, o):
            return nn.Sequential(nn.Conv2d(i, o, 3, padding=1), nn.BatchNorm2d(o), nn.ReLU(),
                                 nn.Conv2d(o, o, 3, padding=1), nn.BatchNorm2d(o), nn.ReLU())
        self.down = nn.ModuleList([block(1, c[0]), block(c[0], c[1]), block(c[1], c[2]), block(c[2], c[3])])
        self.up = nn.ModuleList([block(c[3] + c[2], c[2]), block(c[2] + c[1], c[1]), block(c[1] + c[0], c[0])])
        self.out = nn.Conv2d(c[0], 1, 1)

    def forward(self, x):
        skips = []
        for k, d in enumerate(self.down):
            x = d(x)
            if k < len(self.down) - 1:
                skips.append(x)
                x = F.max_pool2d(x, 2)
        for u in self.up:
            s = skips.pop()
            x = F.interpolate(x, size=s.shape[2:], mode='bilinear', align_corners=False)
            x = u(torch.cat([x, s], dim=1))
        return self.out(x)


def pick_peaks(prob, tissue, n=N_POINTS, min_dist=MIN_DIST):
    """Greedy peaks of the landmark map, at least min_dist apart, inside the tissue."""
    p = np.where(tissue, prob, -np.inf)
    chosen = []
    for k in np.argsort(p, axis=None)[::-1]:
        if not np.isfinite(p.flat[k]) or len(chosen) >= n:
            break
        y, x = divmod(int(k), p.shape[1])
        if all((y - cy) ** 2 + (x - cx) ** 2 >= min_dist ** 2 for cy, cx in chosen):
            chosen.append((y, x))
    return np.array(chosen, dtype=float).reshape(-1, 2)


@torch.no_grad()
def landmarks(net, plane_u8):
    img = plane_u8.astype(np.float32) / 255
    prob = torch.sigmoid(net(torch.tensor(img, device=device())[None, None]))[0, 0].cpu().numpy()
    return pick_peaks(prob, plane_u8 > TISSUE)


# ------------------------------------------------------------------- matcher

class MatchNet(nn.Module):
    """Small ConvNet: (atlas patch, registered-section patch) -> shift (dy, dx) in px."""

    def __init__(self):
        super().__init__()
        c = [2, 32, 64, 96, 128, 160]
        layers = []
        for i in range(5):
            layers += [nn.Conv2d(c[i], c[i + 1], 3, padding=1), nn.BatchNorm2d(c[i + 1]), nn.ReLU(),
                       nn.Conv2d(c[i + 1], c[i + 1], 3, padding=1), nn.BatchNorm2d(c[i + 1]), nn.ReLU(),
                       nn.MaxPool2d(2)]
        self.features = nn.Sequential(*layers)                          # 96 -> 3
        self.head = nn.Sequential(nn.Flatten(), nn.Linear(160 * 9, 128), nn.ReLU(), nn.Linear(128, 2))

    def forward(self, x):
        return self.head(self.features(x))


def _patches(warped, atlas01, pts):
    half = STORE // 2
    Wp = np.pad(warped, half); Ap = np.pad(atlas01, half)
    yx = np.rint(pts).astype(int)
    atl = np.stack([(Ap[y:y + STORE, x:x + STORE] * 255).astype(np.uint8) for y, x in yx])
    hist = np.stack([(np.clip(Wp[y:y + STORE, x:x + STORE], 0, 1) * 255).astype(np.uint8) for y, x in yx])
    return atl, hist


@torch.no_grad()
def correct(net, warped, atlas01, pts):
    """The learned shift of every landmark (atlas frame, relative to the landmark pixel), and its spread."""
    atl, hist = _patches(warped, atlas01, pts)
    c = SHIFT
    a = atl[:, c:c + P, c:c + P]
    answers = []
    for t in TTA:
        h = hist[:, c + t[0]:c + t[0] + P, c + t[1]:c + t[1] + P]
        x = torch.tensor(np.stack([a, h], axis=1).astype(np.float32) / 255, device=device())
        answers.append(net(x).cpu().numpy() + np.array(t))
    answers = np.stack(answers)                                          # 9, n, 2
    d = answers.mean(axis=0)
    spread = np.sqrt(((answers - d) ** 2).sum(axis=2).mean(axis=0))
    return d, spread


# ------------------------------------------------------------------ the chain

def load_models(landmark_path, matcher_path):
    dev = device()
    lnet, mnet = LandmarkNet().to(dev), MatchNet().to(dev)
    lnet.load_state_dict(torch.load(landmark_path, map_location=dev, weights_only=False)['state'])
    mnet.load_state_dict(torch.load(matcher_path, map_location=dev, weights_only=False)['state'])
    return lnet.eval(), mnet.eval()


def propose(sections_rgb, planes_u8, lnet, mnet):
    """Control points for every section.

    sections_rgb  list of (H, W, 3) uint8 volume_for_inspection pages
    planes_u8     list of (H, W) uint8 atlas planes, one per section, as the GUI draws them
    returns one dict per section: atlas (n, 2), hist (n, 2), spread (n,), low (n,) bool
    """
    secs = [section_image(s) for s in sections_rgb]
    atl01 = [p.astype(np.float32) / 255 for p in planes_u8]
    theta, field = register(secs, atl01)
    out = []
    for i, (sec, plane) in enumerate(zip(secs, planes_u8)):
        h, w = plane.shape
        a = landmarks(lnet, plane)
        if len(a) == 0:
            out.append(dict(atlas=np.zeros((0, 2)), hist=np.zeros((0, 2)), spread=np.zeros(0), low=np.zeros(0, bool)))
            continue
        warped = warp_to_atlas(sec, theta[i], field[i])
        d, spread = correct(mnet, warped, atl01[i], a)
        p = map_points(theta[i], field[i], np.rint(a) + d, h, w)
        ok = (p[:, 0] >= 0) & (p[:, 0] < h) & (p[:, 1] >= 0) & (p[:, 1] < w)
        a, p, spread = a[ok], p[ok], spread[ok]
        low = spread > np.quantile(spread, 1 - LOW_CONFIDENCE) if len(spread) else np.zeros(0, bool)
        out.append(dict(atlas=a, hist=p, spread=spread, low=low))
    return out
