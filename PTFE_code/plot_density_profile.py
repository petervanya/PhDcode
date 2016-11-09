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
plt.axis("off")
outname = fname.rstrip("out") + "png"
plt.savefig(outname)
print("2d density profile plot saved in %s." % outname)


