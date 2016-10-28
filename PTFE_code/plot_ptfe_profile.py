#!/usr/bin/env python
"""Usage:
    plot_ptfe_profile.py <fileW> <fileB> [--boxsize <L> --pdf]

[AD HOC] script to format plotting of normalised density profiles
for water and PTFE backbone. (Previously also plotted sulfonic acid groups.)

Options:
    --boxsize <L>       Box size [default: 20]
    --pdf               Plot as pdf

25/02/16
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys, os
from docopt import docopt

rc = 8.14e-10
matplotlib.rcParams.update({'font.size': 16})

args = docopt(__doc__)
fileW = args["<fileW>"]
fileB = args["<fileB>"]
if not os.path.isfile(fileW) or not os.path.isfile(fileB):
    print("One or two files not found.")
    sys.exit()
    
L = float(args["--boxsize"])
lw = 3.0

A = np.loadtxt(fileB)
A[:, 0] *= rc*1e9   # convert to nm
A[:, 1] /= sum(A[:, 1])
plt.plot(A[:, 0], A[:, 1], "red", label="backbone", linewidth=lw)

A = np.loadtxt(fileW)
A[:, 0] *= rc*1e9
A[:, 1] /= sum(A[:, 1])
plt.plot(A[:, 0], A[:, 1], "blue", label="water", linewidth=lw)

plt.xlabel("$x$ (nm)", fontsize=20)
plt.xlim([0.5, np.max(A[:, 0]) - 0.5])
plt.legend(loc="best")

imgname = "profile_1d.png"
if args["--pdf"]: imgname = "profile_1d.pdf"
plt.savefig(imgname)
print("Plot saved in", imgname)


