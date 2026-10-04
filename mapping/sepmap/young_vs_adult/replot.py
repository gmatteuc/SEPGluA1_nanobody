"""Redraw the young against adult slice figures from the maps run_compare.py saved.

A change to the figures does not need the comparison again: the maps come from
volumes_ccf20.npz and the figures from young_vs_adult.compare.draw_figures, the
same code as in a full run.

Run by run_replot.py.
"""

import numpy as np

from sepmap.young_vs_adult.compare import OUT, draw_figures


def main() -> None:
    """Redraw the slice figures from the volumes young_vs_adult.compare saved."""
    z = np.load(OUT / "volumes_ccf20.npz")
    draw_figures({k: z[k] for k in z.files}, OUT)
    print("figures redrawn in", OUT)
