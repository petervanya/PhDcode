#!/usr/bin/env python
"""Usage:
    plot_ptfe_profile.py <fileW> <fileB> [--boxsize <L> --slab <sl> --pdf]

[AD HOC] Plot normalised 1d density profiles of water and PTFE backbone.
Can add lines denoting slab width.

Options:
    --boxsize <L>       Box size in DPD units [default: 40]
    --slab <sl>         Slab width in nm
    --pdf               Plot as pdf

25/02/16
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys, os
from docopt import docopt

r_DPD = 8.14e-10
matplotlib.rcParams.update({'font.size': 16})

args = docopt(__doc__)
fileW, fileB = args["<fileW>"], args["<fileB>"]
if not os.path.isfile(fileW) or not os.path.isfile(fileB):
    sys.exit("One or two files not found. Input: %s, %s" % (fileW, fileB))

print("Plotting: %s | %s" % (fileW, fileB))
    
L = float(args["--boxsize"]) * r_DPD * 1e9

if args["--slab"]:
    D = float(args["--slab"])
    plt.axvspan(0.0, (L - D) / 2, facecolor='black', alpha=0.1)
    plt.axvspan((L + D) / 2, L, facecolor='black', alpha=0.1)

A = np.loadtxt(fileB)
A[:, 0] *= r_DPD * 1e9   # convert to nm
A[:, 1] /= sum(A[:, 1])
plt.plot(A[:, 0], A[:, 1], "red", label="backbone", linewidth=2)

A = np.loadtxt(fileW)
A[:, 0] *= r_DPD * 1e9
A[:, 1] /= sum(A[:, 1])
plt.plot(A[:, 0], A[:, 1], "blue", label="water", linewidth=4)

plt.xlabel("$x$ (nm)", fontsize=20)
plt.xlim([0.5, np.max(A[:, 0]) - 0.5])

imgname = "profile_1d.pdf" if args["--pdf"] else "profile_1d.png"
plt.savefig(imgname)
print("Plot saved in", imgname)


