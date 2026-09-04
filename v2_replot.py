"""
Redraw the v2 slice figures from the saved volumes, without redoing the
ten-minute transform. The figure code lives in v2_compare.draw_figures.

  D:\\sep_histology\\code\\tools\\venv_atlas\\Scripts\\python.exe v2_replot.py
"""

import os

import numpy as np

from v2_compare import OUT, draw_figures

if __name__ == '__main__':
    z = np.load(os.path.join(OUT, 'volumes_ccf20.npz'))
    draw_figures({k: z[k] for k in z.files}, OUT)
    print('figures redrawn in', OUT)
