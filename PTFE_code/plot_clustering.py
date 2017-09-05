#!/usr/bin/env python
"""Usage: plot_clustering.py <fname> [--rc <rc> --d <d> --l <l> --bulk --pdf]

Plot the clustering map from the density file
with title giving slab width and water uptake or "bulk".

Options:
    --d <d>        Slab width [default: 5]
    --l <l>        Water uptake [default: 4]
    --rc <rc>      Cutoff [default: 0.3]

18/05/17
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys
from docopt import docopt

matplotlib.rcParams.update({'font.size': 24})


args = docopt(__doc__)
d = int(args["--d"])
l = int(args["--l"])
rc = float(args["--rc"])
fname = args["<fname>"]
try:
    A = np.loadtxt(fname)
    N = len(A)
    B = np.zeros((N, N)).astype(int)
    B[A > rc] = 1
except FileNotFoundError:
    sys.exit("File %s not found." % fname)

ext = ".png"
if args["--pdf"]:
    ext = ".pdf"
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

fig = plt.figure()
ax = plt.axes()
plt.spy(B)
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
plt.title("$d = %i \,\,\\rm{nm}, \quad \lambda = %i$" % (d, l))
if args["--bulk"]:
    plt.title("Bulk, $\quad \lambda = %i$" % l)

figname = "clustering_d%i_l%i" % (d, l) + ext
if args["--bulk"]:
    figname = "clustering_l%i" % (l,) + ext
plt.savefig(figname, bbox_inches="tight")
print("Plot saved in %s." % figname)


