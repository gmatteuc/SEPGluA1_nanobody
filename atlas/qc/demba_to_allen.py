"""Carry the DeMBA P20 annotation into adult CCF space, for check_demba_to_allen.

The young brains register to the DeMBA P20 template and the adults to the CCF, so
voxelwise work across the two needs a way to carry P20 data into adult space. CCF
Translator (brainglobe_ccf_translator) provides it, built from the deformations
that made DeMBA (Carey 2025): the annotation goes from P20 to P56 and from the
demba_dev_mouse space to allen_mouse, as a segmentation, so by nearest neighbour,
since the mean of two region ids is not a region. DeMBA's labels are the CCF's
warped onto the P20 template, so the result should recover the adult annotation
almost exactly; check_demba_to_allen measures how well, against the crop alignment
alone.

The two spaces order their axes differently (the package's metadata):

    demba_dev_mouse   (ML, DV, AP), 11400 x 8000 x 14100 um
    allen_mouse       (AP, DV, ML), 13200 x 8000 x 11400 um

The atlas folder stores the annotation as (AP, DV, ML), so it is permuted before
it goes in; it comes back in the allen order, on the CCF grid at 20 um.

Writes annotation_in_allen_space_20um.nii.gz in atlas_demba_p20 under the data
root, which check_demba_to_allen reads. Runs in tools\\venv_atlas, from the code
root; the deformation fields are downloaded on first use.

    python atlas\\qc\\demba_to_allen.py
"""

import argparse
import importlib.util
from pathlib import Path

import nibabel as nib
import numpy as np
from brainglobe_ccf_translator import Volume


def _build_script():
    """build_demba_atlas.py, one folder up, loaded for its data root and guards."""
    path = Path(__file__).resolve().parents[1] / "build_demba_atlas.py"
    spec = importlib.util.spec_from_file_location("build_demba_atlas", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# the data root by the rule and guards of get_paths.m, as build_demba_atlas.py
# (outside the sepmap package too) works them out, taken from it, not repeated
DATA = _build_script().DATA
DEMBA_ANN = DATA / "atlas_demba_p20" / "annotation_10.nii.gz"
CCF_ANN = DATA / "atlas" / "annotation_10.nii.gz"
OUT = DATA / "atlas_demba_p20" / "annotation_in_allen_space_20um.nii.gz"

# the DeMBA age and the voxel size of its volumes (named _10 for LightSuite, 20 um)
AGE_PND = 20
VOXEL_UM = 20


def load_demba_annotation() -> np.ndarray:
    """The DeMBA P20 annotation as the atlas folder stores it, (AP, DV, ML)."""
    ann = np.asarray(nib.load(str(DEMBA_ANN)).dataobj)
    if ann.shape != (705, 400, 570):
        raise ValueError(
            f"{DEMBA_ANN} is {ann.shape}, expected (705, 400, 570): DeMBA P20 at "
            "20 um, (AP, DV, ML), as build_demba_atlas.py writes it"
        )
    return ann


def to_allen_space(ann: np.ndarray) -> np.ndarray:
    """The annotation carried from DeMBA P20 to the CCF at P56, in (AP, DV, ML)."""
    vol = Volume(
        values=np.transpose(ann, (2, 1, 0)),
        space="demba_dev_mouse",
        voxel_size_micron=VOXEL_UM,
        age_PND=AGE_PND,
        segmentation_file=True,
    )
    vol.transform(target_age=56, target_space="allen_mouse")
    return np.asarray(vol.values)


def check_ccf_grid(out: np.ndarray) -> None:
    """Stop unless the result is on the CCF grid at 20 um, as the check compares it."""
    ccf_shape = nib.load(str(CCF_ANN)).shape
    expected = tuple((n + 1) // 2 for n in ccf_shape)
    if out.shape != expected:
        raise ValueError(
            f"the carried annotation is {out.shape}, but the CCF at 20 um is "
            f"{expected}; check_demba_to_allen needs the two on one grid"
        )


def main():
    """Carry the annotation, check its grid, and write it in mm units."""
    # the annotation, permuted into the demba_dev_mouse order and carried
    ann = load_demba_annotation()
    print(f"DeMBA P20 annotation {ann.shape}, {len(np.unique(ann))} labels", flush=True)
    print("carrying P20 -> P56 and demba_dev_mouse -> allen_mouse ...", flush=True)
    out = to_allen_space(ann)

    # on the CCF grid at 20 um, then written as the atlas folders are
    check_ccf_grid(out)
    img = nib.Nifti1Image(out.astype(np.uint16), np.diag([0.02, 0.02, 0.02, 1.0]))
    img.header.set_xyzt_units("mm")
    nib.save(img, str(OUT))
    print(f"wrote {OUT} {out.shape}, {len(np.unique(out))} labels")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="the DeMBA P20 annotation carried into adult CCF space"
    )
    parser.parse_args()
    main()
