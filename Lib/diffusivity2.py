#!/usr/bin/env python
"""Usage:
    diffusivity2.py 1d <frames> [--bt <bt> --nb <nb> --si --mols-bead <mb>]
    diffusivity2.py 2d <frames> [--bt <bt> --nb <nb> --si --mols-bead <mb>]
    diffusivity2.py 3d <frames> [--bt <bt> --nb <nb> --si --mols-bead <mb>]

Compute diffusivity (Einstein coeff) 
in bulk and slab materials for a given bead type.
  D = <r^2> / (n t), n = 2 (1d), 4 (2d), 6 (3d)
Averating through at least 50 frames recommended (navg).
Options:
    --bt <bt>         Bead type [default: 1]
    --nb <nb>         Number of beads to average over [default: 1]
    --si              Convert to SI units
    --mols-bead <mb>  Rescale by number of molecules in a bead [default: 6]

pv278@cam.ac.uk, 27/10/16
"""
import numpy as np
from math import sqrt
import time
import os, sys, glob
import matplotlib.pyplot as plt
from docopt import docopt
from lmp_lib import read_xyzfile

AMU = 1.66e-27
m0 = 6 * 18 * AMU
kB = 1.381e-23
T = 300.0
r_DPD = 8.14e-10
tau_DPD = sqrt(m0 * r_DPD**2 / (kB * T))


def diff_xyz(xyz, L, cut=5.0):
    """Create matrix of differences"""
    N, dim = xyz.shape
    d_xyz = np.zeros((N-1, dim))
    for i in range(N-1):
        d_xyz[i] = xyz[i+1] - xyz[i]

    d_xyz[d_xyz > cut] -= L
    d_xyz[d_xyz < -cut] += L
    return d_xyz


def gen_rsq(xyz, L, cut=5.0):
    """Plot <r^2> vs time averaged over many particles"""
    N = len(xyz)
    rsq = np.zeros(N)
    d_xyz = diff_xyz(xyz, L, cut=cut)
    for i in range(1, N):
        rsq[i] = np.sum(np.sum(d_xyz[0:i], 0)**2)
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
    if bt not in bts:
        sys.exit("Bead type %i not found in frame." % bt)
    return


def parse_control():
    """Extract timestep from CONTROL file"""
    if not os.path.isfile("CONTROL"):
        sys.exit("CONTROL file not found.")
    f = open("CONTROL", "r").readlines()
    for line in f:
        if "volume" in line:
            L = float(line.split()[1])**(1. / 3)
        if "timestep" in line:
            dt = float(line.split()[1])
        if "trajectory" in line:
            freq = int(line.split()[2])
    return dt, freq, L


if __name__ == "__main__":
    args = docopt(__doc__)
    frames = glob.glob(args["<frames>"])
    frames.sort()
    Nf = len(frames)
    bt = int(args["--bt"])
    check_beadtypes(bt, frames[0])
    dt, freq, L = parse_control()
    tau = dt * freq
    Nb = int(args["--nb"])
    cut = 5.0  # THINK THROUGH
    mols_bead = int(args["--mols-bead"])

    print("===== Diffusivity =====")
    print("Frames: %i | dt: %.3f | steps/frame: %i" % (Nf, dt, freq))
    print("Bead type: %i | Number of beads: %i | Molecules/bead: %i" % \
            (bt, Nb, mols_bead))
    print("PBC cutoff: %.2f" % cut)

    # Reading frames
    xyzs = np.zeros((Nb, 3, Nf))   # (beads, 3, frames)
    for i in range(Nf):
        xyzs[:, :, i] = read_from_top(frames[i], Nb, btype=bt)
    Rsq = np.zeros(Nf)
    t = np.arange(0, Nf * tau, tau)

    if args["1d"]:
        print("1d diffusivity")
        all_Rs = np.zeros((Nf, 3))
        Ds = np.zeros(3)
        for ax in range(3):
            print("Axis: %i. " % ax, end="")

            ti = time.time()
            for i in range(Nb):
                Rsq = gen_rsq(xyzs[i, (ax,), :].T, L)
                all_Rs[:, ax] += Rsq / Nb
            tf = time.time()

            print("Time: %.2f s." % (tf - ti))

        for j in range(3):
            Ds[j] = np.polyfit(t, all_Rs[:, j], 1)[0] / 2.0 * mols_bead
        print("Diffusivities (x, y, z):", Ds)
        if args["--si"]:
            print("In SI:", Ds * r_DPD**2 / tau_DPD)
        
        outname = "diffusivity_1d_b%i.out" % bt
        hdr = "t, rsq(x), rsq(y), rsq(z)"
        np.savetxt(outname, np.hstack((np.matrix(t).T, all_Rs)), header=hdr)

    if args["2d"]:
        print("2d diffusivity")
        all_Rs = np.zeros((Nf, 3))
        Ds = np.zeros(3)
        planes = [(0, 1), (1, 2), (0, 2)]

        for j in range(3):
            print("Plane: (%i, %i). " % (planes[j][0], planes[j][1]), end="")
            ti = time.time()
            for i in range(Nb):
                Rsq = gen_rsq(xyzs[i, planes[j], :].T, L)
                all_Rs[:, j] += Rsq / Nb
            tf = time.time()
            print("Time: %.2f s." % (tf - ti))

        for j in range(3):
            Ds[j] = np.polyfit(t, all_Rs[:, j], 1)[0] / 4.0 * mols_bead
        print("Diffusivities (xy, yz, xz):", Ds)
        if args["--si"]:
            print("In SI:", Ds * r_DPD**2 / tau_DPD)
        
        outname = "diffusivity_2d_b%i.out" % bt
        hdr = "t, rsq(xy), rsq(yz), rsq(xz)"
        np.savetxt(outname, np.hstack((np.matrix(t).T, all_Rs)), header=hdr)

    if args["3d"]:
        print("3d diffusivity")
        all_Rs = np.zeros(Nf)

        ti = time.time()
        for i in range(Nb):
            Rsq = gen_rsq(xyzs[i, :, :].T, L)
            all_Rs += Rsq / Nb
        tf = time.time()
        print("Time: %.2f s." % (tf - ti))

        D = np.polyfit(t, all_Rs, 1)[0] / 6.0 * mols_bead
        print("Diffusivity:", D)
        if args["--si"]:
            print("In SI:", D * r_DPD**2 / tau_DPD)

        outname = "diffusivity_3d_b%i.out" % bt
        np.savetxt(outname, np.vstack((t, all_Rs)).T)

    print("Output saved in %s." % outname)


