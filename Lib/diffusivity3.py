#!/usr/bin/env python
"""Usage:
    diffusivity3.py <frames> [--bt <bt> --nb <nb> --si --mb <mb> --cut <rc>]
                             [--start <st>]

Compute diffusivity (Einstein coeff) 
in bulk and slab materials for a given bead type.
  D = <r^2> / (n t), n = 2 (1d), 4 (2d), 6 (3d)
Averating through at least 50 frames recommended (navg).

Options:
    --bt <bt>      Bead type [default: 1]
    --nb <nb>      Number of beads for averaging, or 'all' [default: 1]
    --si           Convert to SI units
    --mb <mb>      Rescale by number of molecules in a bead [default: 6]
    --start <st>   Fraction of frames after which to start fitting a line
                   over rsq vs t [default: 0.666]

pv278@cam.ac.uk, 15/11/16
"""
import numpy as np
from math import sqrt
from scipy.optimize import curve_fit
import time, os, sys, glob
from docopt import docopt
from dlms_lib import read_xyzfile

AMU = 1.66e-27
m0 = 6 * 18 * AMU
kB = 1.381e-23
T = 300.0
r_DPD = 8.14e-10
tau_DPD = sqrt(m0 * r_DPD**2 / (kB * T))


def diff_xyz(xyz, L, cut):
    """Create matrix of differences.
    * cut: vector of size 3
    * L: vector of size 3"""
    N, dim = xyz.shape
    d_xyz = np.zeros((N-1, dim))
    for i in range(N-1):
        d_xyz[i] = xyz[i+1] - xyz[i]

#    d_xyz[d_xyz > cut] -= L
#    d_xyz[d_xyz < -cut] += L
    d_xyz = d_xyz - (d_xyz > cut) * L
    d_xyz = d_xyz + (d_xyz < -cut) * L
    return d_xyz


def gen_rsq(xyz, L, cut):
    """Plot <r^2> vs time averaged over many particles"""
    N = len(xyz)
    rsq = np.zeros((N, 3))
    d_xyz = diff_xyz(xyz, L, cut)
    for i in range(1, N):
        rsq[i] = np.sum(d_xyz[0:i], 0)**2
    return rsq


def read_from_top(fname, Nb, btype=1):
    """Read first Nb suplied bead types of from an xyz frame"""
    cnt = 0
    xyz = np.zeros((Nb, 3))
    with open(fname, "r") as f:
        while cnt < Nb:
            line = f.readline()
            if line.split()[0] == str(btype):
                xyz[cnt] = [float(s) for s in line.split()[1:4]]
                cnt += 1
    return xyz


def check_beadtypes(bt, frame):
    """Check if beadtype present in the xyz frame"""
    A = read_xyzfile(frame)
    bts = set(A[:, 0].astype(int))
    Nbf = sum(A[:, 0] == bt)
    if bt not in bts:
        sys.exit("Bead type %i not found in frame." % bt)
    return Nbf


def parse_control():
    """Extract timestep from CONTROL file"""
    if not os.path.isfile("CONTROL"):
        sys.exit("CONTROL file not found.")
    f = open("CONTROL", "r").readlines()
    for line in f:
        if "volume" in line:
#            L = float(line.split()[1])**(1. / 3)
            tmp = line.split()
            if len(tmp) == 2:
                L = float(tmp[-1])**(1/3)
                L = np.array([L, L, L])
            elif len(tmp) == 4:
                L = np.array(list(map(float, tmp[1:4])))
        if "timestep" in line:
            dt = float(line.split()[1])
        if "trajectory" in line:
            freq = int(line.split()[2])
    return dt, freq, L


if __name__ == "__main__":
    args = docopt(__doc__)
    frames = sorted(glob.glob(args["<frames>"]))
    if len(frames) == 0:
        sys.exit("0 frames captured.")
    Nf = len(frames)
    bt = int(args["--bt"])
    dt, freq, L = parse_control()
    tau = dt * freq
    Mb = int(args["--mb"])    # molecules per bead
    Nbf = check_beadtypes(bt, frames[0])
    Nb = args["--nb"]
    if Nb == "all":
        Nb = Nbf
    else:
        Nb = int(Nb)
    if Nb > Nbf:
        sys.exit("Nb larger than number of beads in frame (%i)." % Nbf)

    print("===== Diffusivity =====")
    print("Frames: %i | dt: %.3f | steps/frame: %i" % (Nf, dt, freq))
    print("Bead type: %i | Number of beads: %i | Molecules/bead: %i" % \
            (bt, Nb, Mb))
    print("Box size: %s" % L)

    # Reading frames
    xyzs = np.zeros((Nb, 3, Nf))   # (beads, 3, frames)
    for i in range(Nf):
        xyzs[:, :, i] = read_from_top(frames[i], Nb, btype=bt)
    Rsq = np.zeros((Nf, 3))
    t = np.arange(0, Nf * tau, tau)
    start = float(args["--start"])
    if start > 1.0 or start < 0.0:
        sys.exit("Enter start between 0 and 1.")
    N1 = int(Nf * start)  # start fitting from this point

    # Calculating distances
    ti = time.time()
    for i in range(Nb):
        Rsq += gen_rsq(xyzs[i, :, :].T, L, L/2) / Nb
    tf = time.time()
    print("Time: %.2f s." % (tf - ti))

    # 1D
    outname = "diffusivity_1d_b%i.out" % bt
    hdr = "t, rsq(x), rsq(y), rsq(z)"
    np.savetxt(outname, np.hstack((np.matrix(t).T, Rsq)), \
            fmt="%.6e", header=hdr)

    D = np.zeros(3)
    err = np.zeros(3)
    for i in range(3):
        Cc, Ce = np.polyfit(t[N1:], Rsq[N1:, i], 1, cov=True)
        D[i] = Cc[0] / 2.0 * Mb
        err[i] = np.sqrt(np.diag(Ce))[0] / 2.0 * Mb
    print("\n1d diffusivity  (x,  y,  z):  %.6f  %.6f  %.6f" % tuple(D))
    print("1d fit error    (x,  y,  z):  %.6f  %.6f  %.6f" % tuple(err))
    print("1d <r^2> saved in %s." % outname)
    if args["--si"]:
        print("1d in SI:  %.6e  %.6e  %.6e" % tuple(D * r_DPD**2 / tau_DPD))

    # 2D
    Rsq_2d = np.zeros((Nf, 3))
    Rsq_2d[:, 0] = Rsq[:, 0] + Rsq[:, 1]  # xy
    Rsq_2d[:, 1] = Rsq[:, 1] + Rsq[:, 2]  # yz
    Rsq_2d[:, 2] = Rsq[:, 0] + Rsq[:, 2]  # xz

    outname = "diffusivity_2d_b%i.out" % bt
    hdr = "t, rsq(xy), rsq(yz), rsq(xz)"
    np.savetxt(outname, np.hstack((np.matrix(t).T, Rsq_2d)), \
            fmt="%.6e", header=hdr)

    D = np.zeros(3)
    err = np.zeros(3)
    for i in range(3):
        Cc, Ce = np.polyfit(t[N1:], Rsq_2d[N1:, i], 1, cov=True)
        D[i] = Cc[0] / 4.0 * Mb
        err[i] = np.sqrt(np.diag(Ce))[0] / 4.0 * Mb
    print("\n2d diffusivity (xy, yz, xz):  %.6f  %.6f  %.6f" % tuple(D))
    print("2d fit error   (xy, yz, xz):  %.6f  %.6f  %.6f" % tuple(err))
    print("2d <r^2> saved in %s." % outname)
    if args["--si"]:
        print("2d in SI:  %.6e  %.6e  %.6e" % tuple(D * r_DPD**2 / tau_DPD))

    # 3D
    Rsq_3d = np.sum(Rsq, 1)
    outname = "diffusivity_3d_b%i.out" % bt
    np.savetxt(outname, np.vstack((t, Rsq_3d)).T, fmt="%.6e")

    Cc, Ce = np.polyfit(t[N1:], Rsq_3d[N1:], 1, cov=True)
    D = Cc[0] / 6.0 * Mb
    err = np.sqrt(np.diag(Ce))[0] / 6.0 * Mb
    print("\n3d diffusivity: %.6f" % D)
    print("3d fit error  : %.6f" % err)
    print("3d <r^2> saved in %s." % outname)
    if args["--si"]:
        print("3d in SI:  %.6e" % (D * r_DPD**2 / tau_DPD))


