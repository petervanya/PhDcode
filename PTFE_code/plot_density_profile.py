#!/usr/bin/env python
"""Usage: plot_density_profile.py <fname> [--cmap <cm> --colorbar]

[AD HOC] Plot 2d density profile.
For cmap, check
http://matplotlib.org/examples/color/colormaps_reference.html

Options:
    --cmap <cm>    Color version of the plot [default: gray]
    --colorbar     Add colorbar

08/11/16
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
from docopt import docopt


args = docopt(__doc__)
fname = args["<fname>"]
try:
    A = np.loadtxt(fname)
except FileNotFoundError:
    sys.exit("File %s not found." % fname)
cm = args["--cmap"]

print("Plotting density profile, colormap: %s" % cm)
plt.imshow(A, cmap=cm, vmin=0.0)#, vmax=3.0)
if args["--colorbar"]: plt.colorbar()
ax = plt.axes()
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
figname = fname.rstrip("out") + "png"
plt.savefig(figname, bbox_inches="tight")
print("2d density profile plot saved in %s." % figname)


