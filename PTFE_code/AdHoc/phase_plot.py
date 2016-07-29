#!/usr/bin/env python
"""Usage:
    phase_plot.py <fileB> <fileS> <fileW> [--xmin <xmin> --xmax <xmax>]

Options:
    --xmin <xmin>       Plot starts here [default: 1.13e-8]
    --xmax <xmax>       Plot ends here [default: 1.17e-8]

[AD HOC] Make a figure plot as in Dura presentation.

10/04/16
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import glob, sys
from docopt import docopt


args = docopt(__doc__)
xmin = float(args["--xmin"])
xmax = float(args["--xmax"])

fileB = args["<fileB>"]
fileS = args["<fileS>"]
fileW = args["<fileW>"]

lw = 2
B = np.loadtxt(fileB)
S = np.loadtxt(fileS)
W = np.loadtxt(fileW)
axis = B[:, 0]
norm = B[:, 1] + S[:, 1] + W[:, 1]

plt.fill(axis, B[:, 1]/norm, "r", label="Backbone")
plt.fill(axis, (B[:, 1] + S[:, 1])/norm, "g", label="Sulfonic")
plt.fill(axis, (B[:, 1] + S[:, 1] + W[:, 1])/norm, "b", label="Water")

plotname = "phase_plot.png"
plt.savefig(plotname)
print "Plot saved in", plotname
