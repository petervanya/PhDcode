#!/usr/bin/env python
"""Usage:
    diffusivity.py <frames> [--beads <bn> --dt <dt> --freq <f> --navg <navg>]
                            [--si]

Compute parallel or perpendicular diffusivity (Einstein coeff)
D = <r^2> / (n t), n = 2 in 1d, 4 in 2d
of water in slab Nafion.
TODO
====
* fix which dimension is parr and which perp

ptions:
    --beads <bn>   Number of water beads to consider [default: 1]
    --dt <dt>      Time step [default: 0.03]
    --freq <f>     Now often frames are dumped [default: 100]
    --navg <navg>  Average each point over this many frames [default: 50]
    --si           Convert to SI units

pv278@cam.ac.uk, 18/10/16
"""
import numpy as np
from math import sqrt
import time
import sys, glob
from docopt import docopt
from lmp_lib import read_xyzfile

AMU = 1.66e-27
m0 = 6 * 18 * AMU
kB = 1.381e-23
T = 300.0
r_DPD = 8.14e-10
tau_DPD = sqrt(m0 * r_DPD**2 / (kB * T))


def par_D(Rs, tau, Navg=50):
    """Diffusion in 1d perpendicular to the polymer-electrode interface.
    Input:
    * Rs: matrix of y, z positions (Nf, 2)
    * tau: time step between frames
    * Navg: get one D by averaging over this many frames
    Compute one D by averaging over Navg frames, starting
    at a specific step, then average over the starting steps."""
    Nf = len(Rs)
    tau = freq * dt
    Rs2 = np.sum((Rs - Rs[0])**2, axis=1)
    Ds = np.zeros(Nf - Navg + 1)
    for i in range(Nf - Navg + 1):
        Ds[i] = np.sum(Rs2[i:i+Navg]) / (4 * tau * Navg)
    return np.average(Ds)


def perp_D(Rs, tau, Navg=50):
    """Diffusion in 1d perpendicular to the polymer-electrode interface.
    Input:
    * Rs: vector of x positions of length (Nf)
    * tau: time step between frames
    * Navg: get one D by averaging over this many frames
    Compute one D by averaging over Navg frames, starting
    at a specific step, then average over the starting steps."""
    Nf = len(Rs)
    tau = freq * dt
    Rs2 = (Rs - Rs[0])**2
    Ds = np.zeros(Nf - Navg + 1)
    for i in range(0, Nf - Navg + 1):
        Ds[i] = np.sum(Rs2[i:i+Navg]) / (2 * tau * Navg)
    return np.average(Ds)


def read_from_top(fname, Nb, btype=4):
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


def guess_directions(frame):
    """Guess directions x, y, z (0, 1, 2) from the frame:
    * two parallel to the slab
    * one perpendicular
    """
    A = read_xyzfile(frame)
    xyz = A[A[:, 0] == 4][:, 1:4]
    mins, maxs = np.zeros(3), np.zeros(3)
    diffs = np.zeros(3)
    for i in range(3):
        mins[i], maxs[i] = min(xyz[:, i]), max(xyz[:, i])
    diffs = maxs - mins
    perp = (np.argmin(diffs),)
    parr = tuple(set(np.arange(3)).difference(perp))
    print("Guessing directions. Perp:", perp, "| Parr:", parr)
    return perp, parr


if __name__ == "__main__":
    args = docopt(__doc__)
    frames = glob.glob(args["<frames>"])
    Nf = len(frames)
    if Nf < 50:
        sys.exit("Only %i frames found, at least 50 required." % Nf)
    dt = float(args["--dt"])
    freq = int(args["--freq"])
    tau = dt * freq
    Navg = int(args["--navg"])
    Nb = int(args["--beads"])
    if Nb < 1: 
        sys.exit("Enter number of beads >= 1.")

    print("===== Diffusivity =====")
    print("Frames: %i | dt: %.3f | steps/frame: %i" % (Nf, dt, freq))
    print("Beads to use: %i" % (Nb,))

    perp_dir, parr_dir = guess_directions(frames[-1])

    if Nb == 1:
        Rs = np.zeros((Nf, 3))
        print("Reading files... ", end="")
        ti = time.time()
        for i in range(Nf):
            xyz = read_from_top(frames[i], 1, btype=4) # only W bead
            Rs[i] = xyz
        tf = time.time()
        print("Finished reading. Time: %.2f s." % (tf - ti))

        Do = perp_D(Rs[:, perp_dir], tau, Navg)
        Dp = par_D(Rs[:, parr_dir], tau, Navg)
    else:          # multiple water beads
        Rs = np.zeros((Nb, 3, Nf))   # (beads, 3, frames)
        for i in range(Nf):
            xyz = read_from_top(frames[i], Nb, btype=4)
            Rs[:, :, i] = xyz

        Dps, Dos = np.zeros(Nb), np.zeros(Nb)
        for i in range(Nb):
            Dos[i] = perp_D(Rs[i, perp_dir, :][-1], tau, Navg)
            Dps[i] = par_D(Rs[i, parr_dir, :].T, tau, Navg)
        print("Perp:\n", Dos)
        print("Parr:\n", Dps)
        Dp = np.average(Dps)
        Do = np.average(Dos)

    print("Average over: %i frames | Parallel: %.3e | Perpendicular: %.3e" \
            % (Navg, Dp, Do))
    if args["--si"]:
        Do = Do * r_DPD**2 / tau_DPD
        Dp = Dp * r_DPD**2 / tau_DPD
        print("In SI units (m^2/s):")
        print("Average over: %i frames | Parallel: %.3e | Perpendicular: %.3e" \
                % (Navg, Dp, Do))
             

