#!/usr/bin/env python
"""Usage:
    plot_ptfe_profile.py <fileW> <fileB> [--boxsize <L> --slab <sl> --rc <rc>]
                         [--tex --fmt <fmt> --parse_title]

Plot normalised 1d density profiles of water and PTFE backbone.
Can add lines denoting slab width.
Possible lengthscale:
* 8.14 AA (Wu et al.)
* 7.1 AA (Yamamoto)

Options:
    --boxsize <L>   Box size in DPD units [default: 40]
    --slab <sl>     Slab width in nm
    --parse_title   Set title from information in file name
    --rc <rc>       DPD length scale in AA [default: 8.14]
    --tex           Use tex for axes
    --fmt <fmt>     Plot format [default: png]

25/02/16, updated 17/10/17
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys, os
from docopt import docopt

matplotlib.rcParams.update({'font.size': 24})

args = docopt(__doc__)
r_DPD = float(args["--rc"]) * 1e-10 #8.14e-10
fileW, fileB = args["<fileW>"], args["<fileB>"]
if not os.path.isfile(fileW) or not os.path.isfile(fileB):
    sys.exit("One or two files not found. Input: %s, %s" % (fileW, fileB))

if args["--tex"]:
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
figfmt = args["--fmt"]

print("Plotting: %s | %s" % (fileW, fileB))
fig, ax = plt.subplots(figsize=(8, 6))
L = float(args["--boxsize"]) * r_DPD * 1e9
if args["--slab"]:
    D = float(args["--slab"])
    plt.axvspan(0.0, (L - D) / 2, facecolor='black', alpha=0.1)
    plt.axvspan((L + D) / 2, L, facecolor='black', alpha=0.1)

A = np.loadtxt(fileB)
A[:, 0] *= r_DPD * 1e9   # convert to nm
plt.plot(A[:, 0], A[:, 1], "red", label="backbone", lw=2)

A = np.loadtxt(fileW)
A[:, 0] *= r_DPD * 1e9
plt.plot(A[:, 0], A[:, 1], "blue", label="water", lw=4)

plt.xlabel("$x$ (nm)")
plt.ylabel("Density")
plt.xlim([0.5, L - 0.5])
plt.ylim([0, 3])
plt.yticks([0, 1, 2, 3])

# set title
confname = ""
if args["--parse_title"]:
    raw = fileW.rstrip(".out").split("_")
    print("Parsing filename:", raw)
    d = float(raw[-2].lstrip("d"))
    l = float(raw[-1].lstrip("l"))
    plt.title("$d = %i\,\, \\rm{nm}, \quad\lambda = %i$" % (d, l))
    confname = "_d%i_l%i" % (d, l)

imgname = "profile_1d" + confname + "." + figfmt
plt.savefig(imgname, bbox_inches='tight')
print("Plot saved in", imgname)


