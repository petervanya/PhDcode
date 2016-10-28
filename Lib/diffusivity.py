#!/usr/bin/env python
"""Usage:
    diffusivity.py 1d <axis> <frames> [--beadtype <bt> --nbeads <nb>]
                                      [--dt <dt> --freq <f> --navg <navg> --si]
    diffusivity.py 2d <plane> <frames> [--beadtype <bt> --nbeads <nb>]
                                       [--dt <dt> --freq <f> --navg <navg> --si]
    diffusivity.py 3d <frames> [--beadtype <bt> --nbeads <nb>]
                               [--dt <dt> --freq <f> --navg <navg> --si]

Compute diffusivity (Einstein coeff) 
in bulk and slab materials for a given bead type.
  D = <r^2> / (n t), n = 2 (1D), 4 (2D), 6 (3D)
Averating through at least 50 frames recommended (navg).

Arguments: 
    <axis>           Axis along which to measure 1d D, 'x', 'y', 'z'
    <plane>          Plane in which to measure 1d D, 'xy', 'yz', 'xz'

Options:
    --beadtype <bt>  Bead type [default: 1]
    --nbeads <nb>    Number of beads to average over [default: 1]
    --dt <dt>        Time step [default: 0.03]
    --freq <f>       Now often frames are dumped [default: 100]
    --navg <navg>    Average each point over this many frames [default: 50]
    --si             Convert to SI units

pv278@cam.ac.uk, 27/10/16
"""
import numpy as np
from math import sqrt
import time
import os, sys, glob
from docopt import docopt
from lmp_lib import read_xyzfile

AMU = 1.66e-27
m0 = 6 * 18 * AMU
kB = 1.381e-23
T = 300.0
r_DPD = 8.14e-10
tau_DPD = sqrt(m0 * r_DPD**2 / (kB * T))


def calc_D(Rs, dim, tau, Navg=50):
    """Diffusion. Input:
    * Rs: matrix positions (Nf, dim)
    * dim: dimension, 1, 2, or 3
    * tau: time step between frames
    * Navg: get one D by averaging over this many frames
    Compute one D by averaging over Navg frames, starting
    at a specific step, then average over the starting steps."""
    Nf = len(Rs)
    Ds = np.zeros(Nf - Navg + 1)
    if Rs.shape[1] != dim:
        sys.exit("Matrix of positions does not have correct dimension.")
    Rs2 = np.sum((Rs - Rs[0])**2, axis=1)
    for i in range(0, Nf - Navg + 1):
        Ds[i] = np.sum(Rs2[i:i+Navg]) / (2 * (dim+1) * tau * Navg)
    return np.average(Ds)


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
        if "timestep" in line:
            dt = float(line.split()[-1])
        if "trajectory" in line:
            freq = int(line.split()[2])
    return dt, freq


if __name__ == "__main__":
    args = docopt(__doc__)
    frames = glob.glob(args["<frames>"])
    bt = int(args["--beadtype"])
    check_beadtypes(bt, frames[0])
    Navg = int(args["--navg"])
    Nf = len(frames)
    if Nf < Navg:
        sys.exit("Only %i frames found, at least %i required." % (Nf, Navg))
#    dt = float(args["--dt"])
#    freq = int(args["--freq"])
    dt, freq = parse_control()
    tau = dt * freq
    Nb = int(args["--nbeads"])
    if Nb < 1:
        sys.exit("Enter number of beads >= 1.")

    print("===== Diffusivity =====")
    print("Frames: %i | dt: %.3f | steps/frame: %i" % (Nf, dt, freq))
    print("Beads to use: %i" % (Nb,))

    # Reading frames
    Rs = np.zeros((Nb, 3, Nf))   # (beads, 3, frames)
    for i in range(Nf):
        xyz = read_from_top(frames[i], Nb, btype=bt)
        Rs[:, :, i] = xyz
    Ds = np.zeros(Nb)

    if args["1d"]:
        dim = 1
        axes = {"x": 0, "y": 1, "z": 2}
        ax = args["<axis>"]
        if ax not in axes.keys():
            sys.exit("Enter 'x', 'y', 'z' for axis.")
        ax = axes[ax]
        ax = tuple((ax,))   # matrix slice through tuple returns a np.matrix
        print("1d diffusivity | axis: %i" % ax)
        for i in range(Nb):
#            print(Rs[i, ax, :].T)
            Ds[i] = calc_D(Rs[i, ax, :].T, dim, tau, Navg)

    if args["2d"]:
        dim = 2
        planes = {"xy": (0, 1), "yz": (1, 2), "xz": (0, 2)}
        pl = args["<plane>"]
        if pl not in planes.keys():
            sys.exit("Enter 'xy', 'yz', 'xz' for plane.")
        pl = planes[pl]
        print("2d diffusivity | plane: (%i, %i)" % (pl[0], pl[1]))
        for i in range(Nb):
#            print(Rs[i, pl, :].T)
            Ds[i] = calc_D(Rs[i, pl, :].T, dim, tau, Navg)

    if args["3d"]:
        dim = 3
        print("3d bulk diffusivity")
        for i in range(Nb):
#            print(Rs[i, :, :].T)
            Ds[i] = calc_D(Rs[i, :, :].T, dim, tau, Navg)
    print(Ds)
    D = np.average(Ds)
    stdD = sqrt(np.average(Ds**2) - D**2)

    print("Average over: %i frames | D: %.3f | std(D): %.3f" \
            % (Navg, D, stdD))
    if args["--si"]:
        D = D * r_DPD**2 / tau_DPD
        stdD = stdD * r_DPD**2 / tau_DPD
        print("In SI units (m^2/s):")
        print("Average over: %i frames | D: %.3e | std(D): %.3e" \
                % (Navg, D, stdD))
             

