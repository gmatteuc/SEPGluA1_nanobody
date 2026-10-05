"""The two hemispheres of a cohort volume folded onto one.

The young against adult maps, videos and close-ups average the left and the
mirrored right hemisphere, which doubles the brains behind a voxel where both
sides were cut and leaves a voxel with one side where only one was. A count map
folds to the larger count of the two sides, the most brains behind either.

Only numpy is imported, so the flatmap environment (young_vs_adult.closeup) can
import this module too. Used by volumes.cohort, which folds each brain before the
cohort statistics behind the videos, and by young_vs_adult.compare, video,
video_compare and closeup; it sits at the top of the package because both
sub-packages use it.
"""

import numpy as np


def fold(v: np.ndarray) -> np.ndarray:
    """Average the two hemispheres of an (AP, DV, ML) volume, ignoring NaN.

    The right half is mirrored onto the left; the result is (AP, DV, ML/2).
    """
    ml = v.shape[2]
    h = ml // 2
    left, right = v[:, :, :h], v[:, :, ml - h :][:, :, ::-1]
    return np.nanmean(np.stack([left, right]), axis=0)


def fold_count(n: np.ndarray) -> np.ndarray:
    """Fold a count map as fold does, keeping the larger count of the two sides."""
    ml = n.shape[2]
    h = ml // 2
    return np.maximum(n[:, :, :h], n[:, :, ml - h :][:, :, ::-1])
