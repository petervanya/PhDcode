#!/usr/bin/env python
"""Usage:
    phase_plot2.py <fB> <fS> <fW> <fE> [--L <L> --slab <d> --pdf --out <out>]

Options:
    --L <L>        Box size in DPD units [default: 40]
    --slab <d>     Slab width in nm [default: 10]
    --out <out>    Output name, no extension [default: phase_plot]

[AD HOC] Make a figure plot as in Dura presentation.

10/04/16, modified 25/12/16
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import glob, sys
from docopt import docopt

rc = 8.14e-10
matplotlib.rcParams.update({'font.size': 24})


args = docopt(__doc__)
if args["--pdf"]:           # use latex fonts
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')


fileB = args["<fB>"]
fileS = args["<fS>"]
fileW = args["<fW>"]
fileE = args["<fE>"]

B = np.loadtxt(fileB)
S = np.loadtxt(fileS)
W = np.loadtxt(fileW)
E = np.loadtxt(fileE)
axis = B[:, 0] * rc / 1e-9
norm = B[:, 1] + S[:, 1] + W[:, 1] + E[:, 1]

By = B[:, 1] / norm
Sy = S[:, 1] / norm
Wy = W[:, 1] / norm
Ey = E[:, 1] / norm

L = float(args["--L"]) * rc / 1e-9
d = float(args["--slab"])
x1 = (L - d) / 2
x2 = (L + d) / 2
dx = L / len(B)

plt.fill_between(axis, 0, By, facecolor="red")
plt.fill_between(axis, By, By + Sy, facecolor="lime")
plt.fill_between(axis, By + Sy, By + Sy + Wy, facecolor="blue")
plt.fill_between(axis, By + Sy + Wy, By + Sy + Wy + Ey, facecolor="gray")
plt.xlim([dx, L - dx])
plt.ylim([0, 1])
plt.axvline(x=x1, color="black")
plt.axvline(x=x2, color="black")
plt.xlabel("Box coordinate (nm)")

ext = "pdf" if args["--pdf"] else "png"
plotname = args["--out"] + "." + ext
plt.savefig(plotname, bbox_inches='tight')
print("Phase plot saved in", plotname)


