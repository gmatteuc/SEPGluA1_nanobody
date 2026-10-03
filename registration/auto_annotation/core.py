"""The automatic control-point annotation: the algorithms, no file I/O.

The settings are those validated in the sandbox repository
SEPGluA1_autoannotation, and the numbers quoted here are its measurements on the
24 hand-annotated brains, with models that never saw the brain being scored.

For each section, given the atlas plane it sits on:

    1. map       image-only registration of the section onto its atlas plane: an
                 affine, then a smooth deformation on a 12 x 17 control grid, both
                 by normalised gradient fields (edge alignment) on the DAPI
                 channel, with no points involved; it lands each point at about
                 8 px, the level of an affine fitted to a human's own other points
    2. where     the landmark network (a small U-Net) reads the atlas plane alone
                 and paints where a human would click; its 32 strongest peaks, at
                 least 20 px apart, are the landmarks
    3. correct   the matcher (a small ConvNet) compares a 96 px atlas patch with
                 the registered section's patch at each landmark and predicts the
                 shift the registration missed, averaged over 9 shifted crops
                 (8.1 to 7.15 px on held-out brains)
    4. pair      the partner point on the section is the registration applied to
                 the landmark plus its shift

Each pair also gets a confidence: how much the 9 shifted answers disagree, in px.
The least confident quarter of each section is flagged for review. Coordinates
are (y, x), 0-based, in the 400 x 570 frame the GUI shows.

Called by cli.py.
"""

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

# registration, as in the sandbox's register_affine.py, register_deform.py and
# full_auto.py; per scale of the affine: (downsample factor, steps, learning rate)
AFFINE_SCALES = ((4, 150, 3e-2), (2, 100, 1e-2), (1, 60, 5e-3))

# Gaussian sigma of the blur at each downsample factor, in pixels of that scale
BLUR = {4: 1.0, 2: 1.5, 1: 2.0}

# control points of the deformation, rows x columns
GRID = (12, 17)

# per scale of the deformation: (downsample factor, steps, learning rate)
DEFORM_SCALES = ((2, 120, 2e-3), (1, 80, 1e-3))

# stiffness of the deformation, the weight of its bending energy
BENDING = 100.0

# sections registered together
BATCH = 48

# landmarks, as in the sandbox's landmark_net.py: how many per section
N_POINTS = 32

# least distance between two landmarks, in px
MIN_DIST = 20

# atlas grey level (0-255) above which is brain
TISSUE = 12

# matcher, as in the sandbox's learn_match.py: the patch the network sees, in px
P = 96

# margin stored around the patch, in px
SHIFT = 12
STORE = P + 2 * SHIFT

# the shifts (dy, dx) of the 9 crops whose answers are averaged, in px
TTA = [(dy, dx) for dy in (-6, 0, 6) for dx in (-6, 0, 6)]

# share of each section's pairs flagged as the least confident
LOW_CONFIDENCE = 0.25


def device() -> str:
    """The torch device: the CPU, so that two runs propose the same points.

    On the GPU the registration's gradient goes through grid_sample, whose
    backward pass has no deterministic CUDA kernel (torch 2.6 refuses it under
    torch.use_deterministic_algorithms), so two runs on MG914 moved the proposed
    points by a median of 0.07 px and up to 13 px. On the CPU two runs give the
    same points; on the 64-core workstation a brain takes about 40 s, against 25 s
    on its GPU.
    """
    return "cpu"


# ===== Images =====


def section_image(rgb: np.ndarray) -> np.ndarray:
    """The DAPI channel of a volume_for_inspection page (H, W, 3), scaled 0-1.

    The 1st and 99.5th percentiles go to 0 and 1, and values beyond them are
    clipped.
    """
    f = rgb.astype(np.float32)[..., 2]
    lo, hi = np.percentile(f, 1), np.percentile(f, 99.5)
    return np.clip((f - lo) / max(hi - lo, 1e-6), 0, 1).astype(np.float32)


def blur(x: torch.Tensor, sigma: float) -> torch.Tensor:
    """Separable Gaussian blur of a (B, 1, H, W) tensor, `sigma` in px (none at 0)."""
    if sigma <= 0:
        return x

    # the kernel reaches 3 sigma on each side
    r = int(3 * sigma)
    k = torch.exp(
        -(torch.arange(-r, r + 1, device=x.device, dtype=x.dtype) ** 2) / (2 * sigma**2)
    )
    k = k / k.sum()
    x = F.conv2d(x, k.view(1, 1, 1, -1), padding=(0, r))
    return F.conv2d(x, k.view(1, 1, -1, 1), padding=(r, 0))


# ===== Registration =====


def ngf(a: torch.Tensor, b: torch.Tensor, m: torch.Tensor) -> torch.Tensor:
    """Normalised gradient fields: do the edges line up, whatever their contrast.

    Per image of the batch, the mean over the mask `m` of the squared cosine
    between the gradients of `a` and `b`: 1 where every edge lines up.
    """

    def grads(x):
        """Differences along y and along x, on a common (H - 1, W - 1) grid."""
        return x[..., 1:, :-1] - x[..., :-1, :-1], x[..., :-1, 1:] - x[..., :-1, :-1]

    ay, ax = grads(a)
    by, bx = grads(b)
    mm = m[..., :-1, :-1]

    def eps(gy, gx):
        """The edge parameter: half the mean gradient magnitude in the mask, >= 1e-4.

        Computed without gradients, as a constant of the loss.
        """
        with torch.no_grad():
            mag = torch.sqrt(gy * gy + gx * gx + 1e-12)
            e = (
                0.5
                * (mag * mm).sum(dim=(1, 2, 3), keepdim=True)
                / mm.sum(dim=(1, 2, 3), keepdim=True).clamp_min(1)
            )
        return e.clamp_min(1e-4)

    ea, eb = eps(ay, ax), eps(by, bx)
    na = torch.sqrt(ay * ay + ax * ax + ea * ea)
    nb = torch.sqrt(by * by + bx * bx + eb * eb)
    dot = (ay * by + ax * bx) / (na * nb)
    return (dot * dot * mm).sum(dim=(1, 2, 3)) / mm.sum(dim=(1, 2, 3)).clamp_min(1)


def register_affine(moving: torch.Tensor, fixed: torch.Tensor) -> torch.Tensor:
    """Affine theta (B, 2, 3), normalised coordinates, atlas -> section.

    Fitted coarse to fine (AFFINE_SCALES) with Adam, on `moving` (the sections)
    warped onto `fixed` (the atlas planes), by ngf over the atlas's tissue.
    """
    dev = moving.device
    b = moving.shape[0]
    eye = torch.tensor([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], device=dev).expand(b, 2, 3)
    delta = torch.zeros(b, 2, 3, device=dev, requires_grad=True)
    for factor, steps, lr in AFFINE_SCALES:
        mv = F.avg_pool2d(moving, factor) if factor > 1 else moving
        fx = F.avg_pool2d(fixed, factor) if factor > 1 else fixed
        mv, fx = blur(mv, BLUR[factor]), blur(fx, BLUR[factor])
        mask = (fx > 0.05).float()
        opt = torch.optim.Adam([delta], lr=lr)
        for _ in range(steps):
            grid = F.affine_grid(eye + delta, fx.shape, align_corners=False)
            warped = F.grid_sample(mv, grid, align_corners=False, padding_mode="zeros")
            loss = (1 - ngf(warped, fx, mask)).sum()
            opt.zero_grad()
            loss.backward()
            opt.step()
    return (eye + delta).detach()


def bending(field: torch.Tensor) -> torch.Tensor:
    """Bending energy of a (B, 2, rows, cols) displacement grid, one value per image.

    The mean squared second difference along y plus the same along x.
    """
    dyy = field[..., 2:, :] - 2 * field[..., 1:-1, :] + field[..., :-2, :]
    dxx = field[..., :, 2:] - 2 * field[..., :, 1:-1] + field[..., :, :-2]
    return (dyy**2).mean(dim=(1, 2, 3)) + (dxx**2).mean(dim=(1, 2, 3))


def upsample(field: torch.Tensor, h: int, w: int) -> torch.Tensor:
    """The (B, 2, rows, cols) displacement grid resampled to (h, w), bicubic."""
    return F.interpolate(field, size=(h, w), mode="bicubic", align_corners=True)


def register_deform(
    moving: torch.Tensor, fixed: torch.Tensor, theta: torch.Tensor
) -> torch.Tensor:
    """A smooth displacement (B, 2, GRID), normalised units, added to the affine.

    Fitted coarse to fine (DEFORM_SCALES) as register_affine is, with the bending
    energy weighted by BENDING added to the loss.
    """
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
            warped = F.grid_sample(
                mv, base + disp, align_corners=False, padding_mode="zeros"
            )
            loss = ((1 - ngf(warped, fx, mask)) + BENDING * bending(field)).sum()
            opt.zero_grad()
            loss.backward()
            opt.step()
    return field.detach()


def register(
    sections: list[np.ndarray], planes: list[np.ndarray]
) -> tuple[torch.Tensor, torch.Tensor]:
    """Register every section onto its atlas plane, BATCH sections at a time.

    `sections` and `planes` are lists of (H, W) float images, 0-1. Returns the
    affine thetas (N, 2, 3) and the displacement grids (N, 2, rows, cols).
    """
    dev = device()
    thetas, fields = [], []
    for s in range(0, len(sections), BATCH):
        mv = torch.tensor(np.stack(sections[s : s + BATCH]), device=dev)[:, None]
        fx = torch.tensor(np.stack(planes[s : s + BATCH]), device=dev)[:, None]
        th = register_affine(mv, fx)
        thetas.append(th)
        fields.append(register_deform(mv, fx, th))
    return torch.cat(thetas), torch.cat(fields)


def warp_to_atlas(
    section: np.ndarray, theta: torch.Tensor, field: torch.Tensor
) -> np.ndarray:
    """The section resampled into the atlas frame by the registration."""
    h, w = section.shape
    grid = F.affine_grid(theta[None], (1, 1, h, w), align_corners=False) + upsample(
        field[None], h, w
    ).permute(0, 2, 3, 1)
    img = torch.tensor(section, device=theta.device)[None, None]
    return F.grid_sample(img, grid, align_corners=False)[0, 0].cpu().numpy()


def map_points(
    theta: torch.Tensor | np.ndarray,
    field: torch.Tensor,
    pts: np.ndarray,
    h: int,
    w: int,
) -> np.ndarray:
    """Atlas (y, x) points through affine + displacement into the section.

    `pts` is (n, 2) in pixels of the (h, w) frame; the result is too.
    """
    theta = theta.cpu().numpy() if torch.is_tensor(theta) else theta

    # pixel centres to normalised coordinates, -1 to 1 (align_corners=False)
    xn = (2 * pts[:, 1] + 1) / w - 1
    yn = (2 * pts[:, 0] + 1) / h - 1
    xa, ya = theta @ np.vstack([xn, yn, np.ones_like(xn)])

    # the displacement at each point, read off the upsampled grid
    full = upsample(field[None], h, w)
    where = torch.tensor(
        np.column_stack([xn, yn]), dtype=torch.float32, device=field.device
    ).view(1, 1, -1, 2)
    d = F.grid_sample(full, where, align_corners=False)[0, :, 0].cpu().numpy()
    return np.column_stack([((ya + d[1] + 1) * h - 1) / 2, ((xa + d[0] + 1) * w - 1) / 2])


# ===== Landmarks =====


class LandmarkNet(nn.Module):
    """Small U-Net: atlas plane -> where a human would click (logits)."""

    def __init__(self, c: tuple[int, ...] = (16, 32, 64, 96)) -> None:
        """Four levels down and three up, with `c` channels per level."""
        super().__init__()

        def block(i, o):
            """Two 3 x 3 convolutions to `o` channels, each with batch norm and ReLU."""
            return nn.Sequential(
                nn.Conv2d(i, o, 3, padding=1),
                nn.BatchNorm2d(o),
                nn.ReLU(),
                nn.Conv2d(o, o, 3, padding=1),
                nn.BatchNorm2d(o),
                nn.ReLU(),
            )

        self.down = nn.ModuleList(
            [block(1, c[0]), block(c[0], c[1]), block(c[1], c[2]), block(c[2], c[3])]
        )
        self.up = nn.ModuleList(
            [block(c[3] + c[2], c[2]), block(c[2] + c[1], c[1]), block(c[1] + c[0], c[0])]
        )
        self.out = nn.Conv2d(c[0], 1, 1)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """The logits of the landmark map, the size of the input plane."""
        skips = []
        for k, d in enumerate(self.down):
            x = d(x)
            if k < len(self.down) - 1:
                skips.append(x)
                x = F.max_pool2d(x, 2)
        for u in self.up:
            s = skips.pop()
            x = F.interpolate(x, size=s.shape[2:], mode="bilinear", align_corners=False)
            x = u(torch.cat([x, s], dim=1))
        return self.out(x)


def pick_peaks(
    prob: np.ndarray,
    tissue: np.ndarray,
    n: int = N_POINTS,
    min_dist: float = MIN_DIST,
) -> np.ndarray:
    """Greedy peaks of the landmark map, at least `min_dist` apart, inside the tissue.

    At most `n` of them, strongest first, as (n, 2) float (y, x).
    """
    p = np.where(tissue, prob, -np.inf)
    chosen = []
    for k in np.argsort(p, axis=None)[::-1]:
        if not np.isfinite(p.flat[k]) or len(chosen) >= n:
            break
        y, x = divmod(int(k), p.shape[1])
        if all((y - cy) ** 2 + (x - cx) ** 2 >= min_dist**2 for cy, cx in chosen):
            chosen.append((y, x))
    return np.array(chosen, dtype=float).reshape(-1, 2)


@torch.no_grad()
def landmarks(net: nn.Module, plane_u8: np.ndarray) -> np.ndarray:
    """Landmarks of a uint8 atlas plane: the network map's peaks inside the brain."""
    img = plane_u8.astype(np.float32) / 255
    prob = (
        torch.sigmoid(net(torch.tensor(img, device=device())[None, None]))[0, 0]
        .cpu()
        .numpy()
    )
    return pick_peaks(prob, plane_u8 > TISSUE)


# ===== Matcher =====


class MatchNet(nn.Module):
    """Small ConvNet: (atlas patch, registered-section patch) -> shift (dy, dx) in px."""

    def __init__(self) -> None:
        """Five levels of two convolutions and a max pool, then a linear head."""
        super().__init__()
        c = [2, 32, 64, 96, 128, 160]
        layers = []
        for i in range(5):
            layers += [
                nn.Conv2d(c[i], c[i + 1], 3, padding=1),
                nn.BatchNorm2d(c[i + 1]),
                nn.ReLU(),
                nn.Conv2d(c[i + 1], c[i + 1], 3, padding=1),
                nn.BatchNorm2d(c[i + 1]),
                nn.ReLU(),
                nn.MaxPool2d(2),
            ]

        # five poolings take the 96 px patch down to 3 x 3
        self.features = nn.Sequential(*layers)
        self.head = nn.Sequential(
            nn.Flatten(), nn.Linear(160 * 9, 128), nn.ReLU(), nn.Linear(128, 2)
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """The predicted shift (dy, dx) of each pair of patches, (n, 2)."""
        return self.head(self.features(x))


def _patches(warped, atlas01, pts):
    """Atlas and registered-section patches of STORE px around each point, uint8."""
    half = STORE // 2
    Wp = np.pad(warped, half)
    Ap = np.pad(atlas01, half)
    yx = np.rint(pts).astype(int)
    atl = np.stack(
        [(Ap[y : y + STORE, x : x + STORE] * 255).astype(np.uint8) for y, x in yx]
    )
    hist = np.stack(
        [
            (np.clip(Wp[y : y + STORE, x : x + STORE], 0, 1) * 255).astype(np.uint8)
            for y, x in yx
        ]
    )
    return atl, hist


@torch.no_grad()
def correct(
    net: nn.Module, warped: np.ndarray, atlas01: np.ndarray, pts: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """The learned shift of every landmark, and its spread.

    The shift is in the atlas frame, relative to the landmark pixel, (n, 2); the
    spread is how much the 9 shifted crops (TTA) disagree, in px, (n,).
    """
    atl, hist = _patches(warped, atlas01, pts)
    c = SHIFT
    a = atl[:, c : c + P, c : c + P]
    answers = []
    for t in TTA:
        h = hist[:, c + t[0] : c + t[0] + P, c + t[1] : c + t[1] + P]
        x = torch.tensor(
            np.stack([a, h], axis=1).astype(np.float32) / 255, device=device()
        )
        answers.append(net(x).cpu().numpy() + np.array(t))

    # 9 crops x n landmarks x 2
    answers = np.stack(answers)
    d = answers.mean(axis=0)
    spread = np.sqrt(((answers - d) ** 2).sum(axis=2).mean(axis=0))
    return d, spread


# ===== The chain =====


def load_models(landmark_path: str, matcher_path: str) -> tuple[LandmarkNet, MatchNet]:
    """The landmark and matcher networks with their trained weights, in eval mode."""
    dev = device()
    lnet, mnet = LandmarkNet().to(dev), MatchNet().to(dev)
    lnet.load_state_dict(
        torch.load(landmark_path, map_location=dev, weights_only=False)["state"]
    )
    mnet.load_state_dict(
        torch.load(matcher_path, map_location=dev, weights_only=False)["state"]
    )
    return lnet.eval(), mnet.eval()


def propose(
    sections_rgb: list[np.ndarray],
    planes_u8: list[np.ndarray],
    lnet: LandmarkNet,
    mnet: MatchNet,
) -> list[dict[str, np.ndarray]]:
    """Control points for every section.

    `sections_rgb` is a list of (H, W, 3) uint8 volume_for_inspection pages and
    `planes_u8` a list of (H, W) uint8 atlas planes, one per section, as the GUI
    draws them. Returns one dict per section: atlas (n, 2), hist (n, 2), spread
    (n,) and low (n,) bool; a pair whose partner falls outside the section is
    left out.
    """
    secs = [section_image(s) for s in sections_rgb]
    atl01 = [p.astype(np.float32) / 255 for p in planes_u8]
    theta, field = register(secs, atl01)
    out = []
    for i, (sec, plane) in enumerate(zip(secs, planes_u8)):
        h, w = plane.shape
        a = landmarks(lnet, plane)
        if len(a) == 0:
            out.append(
                dict(
                    atlas=np.zeros((0, 2)),
                    hist=np.zeros((0, 2)),
                    spread=np.zeros(0),
                    low=np.zeros(0, bool),
                )
            )
            continue

        # the partner points: the learned shift, then the registration
        warped = warp_to_atlas(sec, theta[i], field[i])
        d, spread = correct(mnet, warped, atl01[i], a)
        p = map_points(theta[i], field[i], np.rint(a) + d, h, w)

        # keep the pairs whose partner falls inside the section, and flag the
        # least confident share of them
        ok = (p[:, 0] >= 0) & (p[:, 0] < h) & (p[:, 1] >= 0) & (p[:, 1] < w)
        a, p, spread = a[ok], p[ok], spread[ok]
        if len(spread):
            low = spread > np.quantile(spread, 1 - LOW_CONFIDENCE)
        else:
            low = np.zeros(0, bool)
        out.append(dict(atlas=a, hist=p, spread=spread, low=low))
    return out
