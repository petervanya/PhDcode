#!/usr/bin/env python
"""Usage: 
    diff_from_msd.py <msd> --fit-start <sta> --fit-stop <sto>

Compute diffusivity from the MSD.

Options:
    --fit-start <sta>   Time point to start fitting from in tau [default: 0]
    --fit-stop <sto>    Time point to fit up to in tau [default: 2]

pv278@cam.ac.uk, 07/09/17
"""
import numpy as np
import sys
from docopt import docopt


args = docopt(__doc__)
fname = args["<msd>"]
try:
    A = np.loadtxt(fname)
except FileNotFoundError:
    sys.exit("File %s not found." % fname)

t, Rsq = A[:, 0], A[:, 1:]
Nf = len(t)
dt = t[1] - t[0]
strt = float(args["--fit-start"])
stop = float(args["--fit-stop"])
if stop < 0.0 or stop > max(t):
    sys.exit("Enter fitting stop between 0 and %.1f." % max(t))
if strt < 0.0 or strt > stop:
    sys.exit("Enter fitting start between 0 and stop %.1f." % stop)

print("===== Fitting MSD =====")
print("dt: %.3f | Frames: %i" % (dt, Nf))

if "_1d" in fname:
    print("Fitting 1d MSDs")
    D, err = np.zeros(3), np.zeros(3)
    for i in range(3):
        Cc, Ce = np.polyfit(t[(t >= strt) & (t <= stop)], 
                Rsq[(t >= strt) & (t <= stop), i], 1, cov=True)
        D[i] = Cc[0] / 2.0
        err[i] = np.sqrt(np.diag(Ce))[0] / 2.0
    print("1d diffusivity  (x,  y,  z):  %.6f  %.6f  %.6f" % tuple(D))
    print("1d fit error    (x,  y,  z):  %.6f  %.6f  %.6f" % tuple(err))

if "_2d" in fname:
    print("Fitting 2d MSDs")
    D, err = np.zeros(3), np.zeros(3)
    for i in range(3):
        Cc, Ce = np.polyfit(t[(t >= strt) & (t <= stop)], 
                Rsq[(t >= strt) & (t <= stop), i], 1, cov=True)
        D[i] = Cc[0] / 4.0
        err[i] = np.sqrt(np.diag(Ce))[0] / 4.0
    print("2d diffusivity (xy, yz, xz):  %.6f  %.6f  %.6f" % tuple(D))
    print("2d fit error   (xy, yz, xz):  %.6f  %.6f  %.6f" % tuple(err))

if "_3d" in fname:
    t, Rsq = A[:, 0], A[:, 1]
    print("Fitting 3d MSD")
    Cc, Ce = np.polyfit(t[(t >= strt) & (t <= stop)], 
            Rsq[(t >= strt) & (t <= stop)], 1, cov=True)
    D = Cc[0] / 6.0
    err = np.sqrt(np.diag(Ce))[0] / 6.0
    print("3d diffusivity: %.6f" % D)
    print("3d fit error  : %.6f" % err)


