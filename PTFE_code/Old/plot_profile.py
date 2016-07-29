#!/usr/bin/env python
"""Usage:
    plot_profile.py <infiles> <xlbl> <ylbl> <ttl> <savefile> [--norm]

[AD HOC] script to format plotting of water profiles
with normalisation

10/02/16
"""
import numpy as np
import matplotlib.pyplot as plt
import glob, sys
from docopt import docopt


args = docopt(__doc__)
infiles = glob.glob(args["<infiles>"])
if len(infiles) == 0:
    print "No infiles captured, aborting."
    sys.exit()
print infiles
#print len(infiles), "infiles to plot."

for infile in infiles:
    A = np.loadtxt(infile)
    if args["--norm"]:
        A[:, 1] /= sum(A[:, 1])
    
    plt.plot(A[:, 0], A[:, 1], label=infile)

plt.xlabel(args["<xlbl>"])
plt.ylabel(args["<ylbl>"])
plt.title(args["<ttl>"])
plt.savefig(args["<savefile>"])
print "Plot saved in", args["<savefile>"]

