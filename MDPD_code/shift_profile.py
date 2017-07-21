#!/usr/bin/env python
"""Usage:
    shift_profile.py <profile> [--where <s> --cut <rhoc> --plot]

Shift the density profile so that the liquid blob is centered at <s>,
by default the centre of the box (s = 0.5)-
Treat everything below a certain cutoff as effectively "zero" density.

Options:
    --where <s>    Where to centre the liquid blob [default: 0.5]
    --cut <rhoc>   Cutoff below which density is zero [default: 0.5]

20/07/17
"""
import numpy as np
from scipy.optimize import curve_fit
import sys
from docopt import docopt


args = docopt(__doc__)
cut = float(args["--cut"])
pr = args["<profile>"]
try:
    A = np.loadtxt(pr)
except FileNotFoundError:
    sys.exit("File %s not found." % pr)

x, rho = A[:, 0], A[:, 1]
rho_new = np.zeros(len(rho))
dx = x[1] - x[0]
N = len(x)
L = N * dx
s = float(args["--where"])
if s > 1.0 or s < 0.0:
    sys.exit("Enter <s> between 0 and 1.")
Nc = int(N * s)

print("===== Shifting the density profile =====")
print("Center at: %.3f" % s)

Nzero = sum(rho < cut)
if Nzero != 0:
    # check if last element is zero
    if rho[-1] <= cut:
        # case with blob in the middle: find indices with nonzero density
#        Old: align the blob on the left side of the box
#        occ = np.argwhere(rho > cut).T[0]
#        start = occ[1]
#        rho_new = np.roll(rho, -(start-1))
        occ = np.argwhere(rho > cut).T[0]
        i_mid = occ[len(occ) // 2]
        rho_new = np.roll(rho, Nc - i_mid)
    else:                
        # case with blob across the boundary
        # index of first element on the right with nonzero density
        emp = np.argwhere(rho <= cut).T[0]
        ifirst = np.max(emp) + 1
        Nr = N - ifirst + 1
        rho_new = np.roll(rho, Nr)      # align the blob on the left
        occ = np.argwhere(rho_new > cut).T[0]
        i_mid = occ[len(occ) // 2]
        rho_new = np.roll(rho_new, Nc - i_mid)
else:
    sys.exit("Homogeneous liquid, no need to shift.")

outname = args["<profile>"].rstrip(".out") + "_shift.out"
np.savetxt(outname, np.c_[x, rho_new])
print("Shifted profile saved in %s." % outname)

if args["--plot"]:
    import matplotlib.pyplot as plt
    plt.plot(x, rho, label="old")
    plt.plot(x, rho_new, label="new")
    plt.legend()
    plt.xlim([0, L])
    plt.ylim(bottom=0)
    plt.show()


