"""Per-brain volumes from the registered stacks, each on the atlas of its own age.

The normalisation of the plasticity chain (run_normalise_groups) would get a
cross-age question wrong, in order of damage:

    - its tissue mask comes from the nano intensity (select_background_pixels),
      which flags a whole plane as background when a section covers only part
      of it, and as tissue when the plane is empty;
    - its per-mouse map is affine (slope + intercept), which does not preserve
      ratios between regions;
    - it takes abs() of the normalised values, and cohort means over voxels
      where some mice have no tissue at all.

This module instead works per mouse, at 20 um (2x2x2 block mean of the
registered 10 um-equivalent grid, over the 10 um voxels a section was imaged at:
a voxel the registration left at 0 is missing, not dark, and enters no mean):

    tissue  auto channel above its off-tissue level (median + 4 MAD), and the
             registered nano non-zero, and inside the atlas brain. The nano
             intensity never enters the mask.
    sig      nano minus the mouse's scalar off-tissue background (raw counts)
    auto     auto minus its own off-tissue background
    sep      the SEP (green) channel, same treatment, where run_add_sep_channel
             has carried it into registered space. It was meant to make nano/sep
             read as surface per unit receptor expressed, against nano/auto's
             surface per unit tissue. It does not: adult.sep_channel_check finds
             this channel dominated by autofluorescence in fixed, mounted tissue.
             Kept because it is a real measurement and the check needs it, but
             nano/sep is not a surface fraction.

Plus two scalars: the isocortex mean of sig (the pure-scale cortex reference of
the cref reading downstream) and the backgrounds.

Off-tissue samples are taken from planes where the atlas says brain and the
registered nano is non-zero on more than half of it (planes a section
reached), outside the atlas brain mask, so a degenerate background mask
cannot pollute them.

Writes comparisons_v2/per_mouse/<mouse>.npz under the data root (sig, auto, sep:
float16, NaN where none of a block's 10 um voxels was imaged; tissue: bool;
scalars). Nothing under comparisons/ is touched.

Run by run_per_mouse.py.
"""

import csv
import time
from pathlib import Path

import nibabel as nib
import numpy as np
import tifffile

from sepmap.config import DATA, SETTINGS, mapping_mice

TISSUE = SETTINGS["tissue"]

OUT = DATA / "comparisons_v2" / "per_mouse"
CSV_MAP = DATA / "atlas" / "parcellation_to_parcellation_term_membership.csv"


def atlas_of(row: dict[str, str]) -> str:
    """The atlas a brain is registered to, as cohort_atlas_key.m decides it.

    The adults are on the Allen CCF, a young brain on the DeMBA atlas of its age.
    """
    if row["group"] == "young":
        return f"demba_p{int(row['age_days'])}"
    return "ccf"


# mouse -> (cohort, atlas, group folder under data\), from the cohort table
# (config.mapping_mice), in the route's order. Every brain is read from its own
# registered tiffs, on the atlas of its own age. MG911 is P16 and MG904 P22: same
# recipe, different grid, so nothing here may assume the P20 crop (atlas_grid).
MICE = {r["name"]: (r["mapping_cohort"], atlas_of(r), r["group"]) for r in mapping_mice()}


def atlas_grid(atlas_key: str) -> tuple[Path, tuple[int, int], None]:
    """(annotation folder, AP crop, None) for an atlas key, read from disk.

    The crop, in 10 um planes counted from 1, is whatever build_demba_atlas.py
    measured for that age and wrote to aplims.txt, the same file get_atlas
    reads in MATLAB, so the Python side cannot drift from the atlas a brain was
    registered against. The third element is always None.
    """
    if atlas_key == "ccf":
        return DATA / "atlas", (180, 1079), None
    d = DATA / ("atlas_" + atlas_key)
    lo, hi = (int(v) for v in (d / "aplims.txt").read_text().split())
    return d, (lo, hi), None


def annotation_20(atlas_key: str) -> np.ndarray:
    """Labels on the 20 um (AP, DV, ML) grid of the block-averaged registered stack.

    CCF ships at 10 um and is cropped then halved; the DeMBA atlases already are
    20 um, so they are only cropped. Either way the result is the grid the
    registered volumes land on after the 2x2x2 block mean, for that mouse's own
    age.
    """
    d, (lo, hi), _ = atlas_grid(atlas_key)
    ann = np.asarray(nib.load(d / "annotation_10.nii.gz").dataobj)

    # CCF: 450 x 400 x 570; DeMBA: 497 planes at P20, 478 at P16
    if atlas_key == "ccf":
        return ann[lo - 1 : hi][::2, ::2, ::2]
    return ann[lo - 1 : hi]


def structure_terms() -> tuple[dict[int, str], dict[int, str], dict[int, str]]:
    """Structure names and acronyms, and division acronyms, by parcellation index.

    Read from the parcellation term membership table (CSV_MAP). The layers of an
    area are separate indices that share one structure name, the key of the
    region tables.
    """
    names, acro, divi = {}, {}, {}
    with open(CSV_MAP, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            idx = int(row["parcellation_index"])
            if row["parcellation_term_set_name"] == "structure":
                names[idx] = row["parcellation_term_name"]
                acro[idx] = row["parcellation_term_acronym"]
            elif row["parcellation_term_set_name"] == "division":
                divi[idx] = row["parcellation_term_acronym"]
    return names, acro, divi


def isocortex_ids() -> set[int]:
    """The parcellation indices of the isocortex division."""
    iso = set()
    with open(CSV_MAP, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            if (
                row["parcellation_term_set_name"] == "division"
                and row["parcellation_term_acronym"] == "Isocortex"
            ):
                iso.add(int(row["parcellation_index"]))
    return iso


class Source:
    """One registered channel of a brain, read page by page.

    A page is the (AP, DV) plane at one ML index. Young and adult alike are read
    from the brain's own registered tiffs, so every channel of a brain comes
    from one registration and SEP joins them on the same footing. The obvious
    alternative for the adults, the cohort-wide nano_4d.mat and auto_4d.mat, is
    harmless for nano (it matches the tiffs voxel for voxel, r = 1.0000) and
    wrong for auto: auto_4d.mat was collected a week before the registration
    the tiffs come from and sits about two 10 um voxels off it (r 0.94-0.99 per
    plane, every adult), so nano over auto would be taken between two
    registrations of the same brain.
    """

    FILES = {
        "nano": ("volume_registered", "chan02_NANO.tiff"),
        "auto": ("volume_registered", "chan03_AUTO.tiff"),
        "sep": ("volume_registered_sep", "chan02_SEP.tiff"),
    }

    def __init__(self, mouse: str, group: str, chan: str) -> None:
        """Open channel `chan` ('nano', 'auto' or 'sep') of `mouse` in `group`."""
        self.path = channel_path(mouse, group, chan)
        self.t = tifffile.TiffFile(self.path)
        self.n = len(self.t.pages)

    def page(self, i: int) -> np.ndarray:
        """Page `i` as float32."""
        return np.asarray(self.t.pages[i].asarray(), dtype=np.float32)


def channel_path(mouse: str, group: str, chan: str) -> Path:
    """The registered tiff of channel `chan` of `mouse`, in group folder `group`."""
    sub, name = Source.FILES[chan]
    return DATA / group / mouse / "lightsuite" / sub / name


def has_sep(mouse: str) -> bool:
    """Whether run_add_sep_channel has carried the brain's SEP channel across."""
    return channel_path(mouse, MICE[mouse][2], "sep").exists()


def block2(a: np.ndarray) -> np.ndarray:
    """2x2 block mean of a 2D page with even dims."""
    return a.reshape(a.shape[0] // 2, 2, a.shape[1] // 2, 2).mean(axis=(1, 3))


def imaged_mean(total: np.ndarray, imaged: np.ndarray) -> np.ndarray:
    """`total` over `imaged`, a block's mean over its imaged voxels; NaN where none."""
    mean = np.full(total.shape, np.nan, np.float32)
    np.divide(total, imaged, out=mean, where=imaged > 0)
    return mean


def block_means(
    mouse: str, group: str, n_ap: int, n_dv: int, n_ml: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray | None, np.ndarray]:
    """The 2x2x2 block mean of nano, auto and SEP, at 20 um, from the registered tiffs.

    Returns (sig, aut, sp, nz): nano, auto and SEP (None without a SEP channel), each
    the mean over the block's imaged 10 um voxels, NaN where none of the 8 was
    imaged, and the fraction of the 8 with nano non-zero. Read two ML pages at a
    time, each 2x2 block-averaged.
    """
    sources = {"nano": Source(mouse, group, "nano"), "auto": Source(mouse, group, "auto")}
    if has_sep(mouse):
        sources["sep"] = Source(mouse, group, "sep")
    means = {name: np.zeros((n_ap, n_dv, n_ml), np.float32) for name in sources}

    # fraction of the 8 fine voxels with nano != 0
    nz = np.zeros((n_ap, n_dv, n_ml), np.float32)

    # 2x2x2 block mean over each channel's imaged 10 um voxels. The registration
    # leaves 0 where no section was imaged; imaged tissue or slide reads about 100
    # counts or more (1st percentile, MG897 and CGF027). So a voxel exactly 0 in a
    # channel is missing there, not dark: the channels agree on these voxels to 1 in
    # 5,000, but the auto channel of MG911 and MG904 is 0 over patches where nano
    # reads 500 to 2000 counts (1.5% and 4% of nano's imaged voxels). A plain mean
    # of the 8 would average them in as zeros, darkening the blocks at a section's
    # edge or a patch's, the off-tissue blocks the backgrounds come from included
    for j in range(n_ml):
        for name, source in sources.items():
            # over the two pages, the sums of the 2x2 block means of the values and
            # of the imaged voxels; their ratio is the mean over the imaged voxels
            total = np.zeros((n_ap, n_dv), np.float32)
            imaged = np.zeros((n_ap, n_dv), np.float32)
            for i in (2 * j, 2 * j + 1):
                v = source.page(i)
                total += block2(v)
                imaged += block2((v != 0).astype(np.float32))
            means[name][:, :, j] = imaged_mean(total, imaged)
            if name == "nano":
                nz[:, :, j] = imaged / 2
    return means["nano"], means["auto"], means.get("sep"), nz


def reached_planes(nz: np.ndarray, brain: np.ndarray, n_ap: int) -> np.ndarray:
    """Planes a section reached: nano non-zero on more than half the atlas brain."""
    min_reached = TISSUE["min_reached"]
    reached = np.array(
        [
            (nz[k][brain[k]] > 0).mean() > min_reached if brain[k].any() else False
            for k in range(n_ap)
        ]
    )
    return reached


def backgrounds(
    sig: np.ndarray,
    aut: np.ndarray,
    nz: np.ndarray,
    brain: np.ndarray,
    reached: np.ndarray,
) -> tuple[np.ndarray, float, float, float]:
    """The off-tissue voxels, and the nano and auto backgrounds and auto MAD there.

    Off tissue is off the atlas brain but imaged, in a plane a section reached;
    each channel's median is over the off-tissue voxels it has a value in. 1.4826
    turns a median absolute deviation into an SD.
    """
    off = (~brain) & (nz > 0) & reached[:, None, None]

    # NaN left out: off tissue is where nano was imaged, and the auto channel was
    # not imaged in some of those blocks (CGF027: 1,862 of 27.0M, MG911: 223,131
    # of 24.5M), where one NaN would make the median NaN and the tissue mask empty
    bg_n = float(np.nanmedian(sig[off]))
    bg_a = float(np.nanmedian(aut[off]))
    mad_a = float(np.nanmedian(np.abs(aut[off] - bg_a))) * 1.4826
    return off, bg_n, bg_a, mad_a


def check_sep_in_tissue(
    mouse: str, group: str, sp: np.ndarray, tissue: np.ndarray
) -> None:
    """Stop if SEP has no value in a tissue voxel, where nano and auto both have one.

    The tissue rule needs nano and auto imaged, but not SEP; the readings smooth
    and warp the channels over the tissue, where one NaN would spread to its
    neighbours (none of the 17 brains of 4 October 2026 has such a voxel).
    """
    n_missing = int(np.isnan(sp[tissue]).sum())
    if n_missing:
        raise ValueError(
            f"{mouse}: SEP was not imaged in {n_missing} tissue voxels at 20 um, where "
            "nano and auto were. Check its registered SEP tiff "
            f"({channel_path(mouse, group, 'sep')}) against run_add_sep_channel.m's "
            "output."
        )


def write_mouse(
    mouse: str,
    sig: np.ndarray,
    aut: np.ndarray,
    tissue: np.ndarray,
    reached: np.ndarray,
    scalars: dict,
    extra: dict,
) -> None:
    """Write <mouse>.npz: the volumes as float16, the mask and the scalars."""
    np.savez_compressed(
        OUT / (mouse + ".npz"),
        sig=sig.astype(np.float16),
        auto=aut.astype(np.float16),
        tissue=tissue,
        reached=reached,
        **scalars,
        **extra,
    )


def main(mice: list[str]) -> None:
    """Write the per-brain file of each of `mice`, with one printed line each."""
    iso = isocortex_ids()
    OUT.mkdir(parents=True, exist_ok=True)
    anns = {}
    for mouse in mice:
        t0 = time.time()
        cohort, atlas_key, group = MICE[mouse]
        if atlas_key not in anns:
            anns[atlas_key] = annotation_20(atlas_key)

        # the atlas brain, (AP, DV, ML) at 20 um, and the block-averaged channels
        ann = anns[atlas_key]
        brain = ann > 0
        n_ap, n_dv, n_ml = ann.shape
        sig, aut, sp, nz = block_means(mouse, group, n_ap, n_dv, n_ml)

        # the planes a section reached, and the backgrounds off tissue
        reached = reached_planes(nz, brain, n_ap)
        off, bg_n, bg_a, mad_a = backgrounds(sig, aut, nz, brain, reached)

        # the tissue mask, on the autofluorescence even with SEP, so that swapping the
        # reference changes only the divisor; without 'reached', which guards only the
        # off-tissue sampling: a plane a section covers at 45% is real data (MG897, MG913)
        tissue = (
            brain & (nz >= TISSUE["min_nonzero"]) & (aut > bg_a + TISSUE["mad_k"] * mad_a)
        )
        sig -= bg_n
        aut -= bg_a

        # the isocortex mean, the scale of the cref reading
        cortex = tissue & np.isin(ann, list(iso))
        cortex_mean = float(sig[cortex].mean())
        lo, hi = atlas_grid(atlas_key)[1]
        extra = {}
        if sp is not None:
            check_sep_in_tissue(mouse, group, sp, tissue)

            # over the off-tissue voxels with SEP imaged, as the other backgrounds
            bg_s = float(np.nanmedian(sp[off]))
            sp -= bg_s
            extra = dict(sep=sp.astype(np.float16), bg_sep=bg_s)
        if extra:
            sep_text = f"bg sep {extra['bg_sep']:5.0f}  "
        else:
            sep_text = "no sep       "
        scalars = dict(
            bg_nano=bg_n,
            bg_auto=bg_a,
            mad_auto=mad_a,
            cortex_mean=cortex_mean,
            cohort=cohort,
            atlas=atlas_key,
            aplims=(lo, hi),
        )
        write_mouse(mouse, sig, aut, tissue, reached, scalars, extra)
        print(
            f"{mouse:20s} {cohort:10s} bg nano {bg_n:6.0f}  "
            f"bg auto {bg_a:5.0f} (mad {mad_a:4.0f})  "
            f"{sep_text}"
            f"tissue {100 * tissue.sum() / brain.sum():5.1f}% of atlas brain, "
            f"planes reached {reached.sum()}/{n_ap}  "
            f"cortex mean {cortex_mean:6.0f}   {time.time() - t0:.0f} s",
            flush=True,
        )
