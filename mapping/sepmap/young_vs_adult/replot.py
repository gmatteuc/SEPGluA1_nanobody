"""
Redraw the v2 slice figures from the saved volumes, without redoing the
ten-minute transform. The figure code lives in young_vs_adult.compare.draw_figures.

  D:\\sep_histology\\code\\tools\\venv_atlas\\Scripts\\python.exe mapping\\run_replot.py
"""

import os

import numpy as np

from sepmap.young_vs_adult.compare import OUT, draw_figures


def main():
    """Redraw the slice figures from the volumes young_vs_adult.compare saved."""
    z = np.load(os.path.join(OUT, 'volumes_ccf20.npz'))
    draw_figures({k: z[k] for k in z.files}, OUT)
    print('figures redrawn in', OUT)
