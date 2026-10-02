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
registered 10 um-equivalent grid):

    tissue  auto channel above its off-tissue level (median + 4 MAD), and the
             registered nano non-zero, and inside the atlas brain. The nano
             intensity never enters the mask.
    sig      nano minus the mouse's scalar off-tissue background (raw counts)
    auto     auto minus its own off-tissue background
    sep      the SEP (green) channel, same treatment, where run_add_sep_channel
             has carried it into registered space. It was meant to make nano/sep
             read as surface per unit receptor expressed, against nano/auto's
             surface per unit tissue. It does not: adult.sep_channel_check finds
             this channel dominated by autofluorescence in fixed, cleared tissue.
             Kept because it is a real measurement and the check needs it, but
             nano/sep is not a surface fraction.

Plus two scalars: the isocortex mean of sig (the pure-scale cortex reference of
the cref reading downstream) and the backgrounds.

Off-tissue samples are taken from planes where the atlas says brain and the
registered nano is non-zero on more than half of it (planes a section
reached), outside the atlas brain mask, so a degenerate background mask
cannot pollute them.

Writes comparisons_v2/per_mouse/<mouse>.npz under the data root (sig, auto, sep:
float16; tissue: bool; scalars). Nothing under comparisons/ is touched.

Run by run_per_mouse.py.
"""

import csv
import time
from pathlib import Path

import nibabel as nib
import numpy as np
import tifffile

from sepmap.config import DATA, SETTINGS

TISSUE = SETTINGS["tissue"]

OUT = DATA / "comparisons_v2" / "per_mouse"
CSV_MAP = DATA / "atlas" / "parcellation_to_parcellation_term_membership.csv"

# mouse -> (cohort, atlas, group folder under data\). Every brain is read from its
# own registered tiffs, on the atlas of its own age. MG911 is P16 and MG904 P22:
# same recipe, different grid, so nothing here may assume the P20 crop (atlas_grid).
MICE = {
    **{
        m: ("young_P20", "demba_p20", "young")
        for m in (
            "MG897_SepGluA_P20",
            "MG903_SepGluA_P20",
            "MG913_SepGluA_P20",
            "MG909_SepGluA_P20",
            "MG910_SepGluA_P20",
        )
    },
    "MG911_SepGluA_P16": ("young_P16", "demba_p16", "young"),
    "MG904_SepGluA_P22": ("young_P22", "demba_p22", "young"),
    **{
        m: ("naive", "ccf", "naive")
        for m in (
            "CGF027_Gria1",
            "CGF028_Gria1",
            "CGF033_Gria1",
            "CGF034_Gria1",
            "CGF035_Gria1",
        )
    },
    **{
        m: ("rws", "ccf", "rws")
        for m in (
            "MG691_Gria1",
            "MG692_Gria1",
            "MG693_Gria1",
            "MG736_Gria1",
            "MG737_Gria1",
        )
    },
}


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
    lo, hi = (int(v) for v in open(d / "aplims.txt").read().split())
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

        # the atlas brain, (AP, DV, ML) at 20 um
        ann = anns[atlas_key]
        brain = ann > 0
        nano, auto = Source(mouse, group, "nano"), Source(mouse, group, "auto")
        sep = Source(mouse, group, "sep") if has_sep(mouse) else None
        n_ap, n_dv, n_ml = ann.shape
        sig = np.zeros((n_ap, n_dv, n_ml), np.float32)
        aut = np.zeros((n_ap, n_dv, n_ml), np.float32)
        sp = np.zeros((n_ap, n_dv, n_ml), np.float32) if sep is not None else None

        # fraction of the 8 fine voxels with nano != 0
        nz = np.zeros((n_ap, n_dv, n_ml), np.float32)

        # 2x2x2 block mean: two ML pages at a time, each 2x2 block-averaged
        for j in range(n_ml):
            s = a = z = g = None
            for i in (2 * j, 2 * j + 1):
                v = nano.page(i)
                u = auto.page(i)
                s = block2(v) if s is None else s + block2(v)
                a = block2(u) if a is None else a + block2(u)
                if z is None:
                    z = block2((v != 0).astype(np.float32))
                else:
                    z = z + block2((v != 0).astype(np.float32))
                if sep is not None:
                    w = sep.page(i)
                    g = block2(w) if g is None else g + block2(w)
            sig[:, :, j] = s / 2
            aut[:, :, j] = a / 2
            nz[:, :, j] = z / 2
            if sep is not None:
                sp[:, :, j] = g / 2

        # planes a section reached: nano non-zero on more than half the atlas brain
        min_reached = TISSUE["min_reached"]
        reached = np.array(
            [
                (nz[k][brain[k]] > 0).mean() > min_reached if brain[k].any() else False
                for k in range(n_ap)
            ]
        )

        # backgrounds from voxels off the atlas brain but imaged; 1.4826 turns a
        # median absolute deviation into an SD
        off = (~brain) & (nz > 0) & reached[:, None, None]
        bg_n = float(np.median(sig[off]))
        bg_a = float(np.median(aut[off]))
        mad_a = float(np.median(np.abs(aut[off] - bg_a))) * 1.4826

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
        if sep is not None:
            bg_s = float(np.median(sp[off]))
            sp -= bg_s
            extra = dict(sep=sp.astype(np.float16), bg_sep=bg_s)
        if extra:
            sep_text = f"bg sep {extra['bg_sep']:5.0f}  "
        else:
            sep_text = "no sep       "
        np.savez_compressed(
            OUT / (mouse + ".npz"),
            sig=sig.astype(np.float16),
            auto=aut.astype(np.float16),
            tissue=tissue,
            reached=reached,
            bg_nano=bg_n,
            bg_auto=bg_a,
            mad_auto=mad_a,
            cortex_mean=cortex_mean,
            cohort=cohort,
            atlas=atlas_key,
            aplims=(lo, hi),
            **extra,
        )
        print(
            f"{mouse:20s} {cohort:10s} bg nano {bg_n:6.0f}  "
            f"bg auto {bg_a:5.0f} (mad {mad_a:4.0f})  "
            f"{sep_text}"
            f"tissue {100 * tissue.sum() / brain.sum():5.1f}% of atlas brain, "
            f"planes reached {reached.sum()}/{n_ap}  "
            f"cortex mean {cortex_mean:6.0f}   {time.time() - t0:.0f} s",
            flush=True,
        )
