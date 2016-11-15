#!/usr/bin/env python
"""Usage: plot_density_profile.py <fname> [--cmap <cm>]

[AD HOC] Plot 2d density profile.
For cmap, check
http://matplotlib.org/examples/color/colormaps_reference.html

Options:
    --cmap <cm>    Color version of the plot [default: gray]

08/11/16
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
from docopt import docopt


args = docopt(__doc__)
fname = args["<fname>"]
A = np.loadtxt(fname)
cm = args["--cmap"]

print("Plotting density profile, colormap: %s" % cm)
plt.imshow(A, cmap=cm, vmin=0.0, vmax=3.0)
#plt.colorbar()
#plt.axis("off")
ax = plt.axes()
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
figname = fname.rstrip("out") + "png"
plt.savefig(figname, bbox_inches="tight")
print("2d density profile plot saved in %s." % figname)


