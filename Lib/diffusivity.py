#!/usr/bin/env python
"""Usage:
    diffusivity.py 1d <ax> <frames> 
                      [--bt <bt> --nbeads <nb> --navg <navg> --si --plot-rsq]
    diffusivity.py 2d <pl> <frames> 
                      [--bt <bt> --nbeads <nb> --navg <navg> --si --plot-rsq]
    diffusivity.py 3d <frames> 
                      [--bt <bt> --nbeads <nb> --navg <navg> --si --plot-rsq]

Compute diffusivity (Einstein coeff) 
in bulk and slab materials for a given bead type.
  D = <r^2> / (n t), n = 2 (1D), 4 (2D), 6 (3D)
Averating through at least 50 frames recommended (navg).

Arguments: 
    <ax>             Axis along which to measure 1d D, 'x', 'y', 'z'
    <pl>             Plane in which to measure 1d D, 'xy', 'yz', 'xz'

Options:
    --bt <bt>        Bead type [default: 1]
    --nbeads <nb>    Number of beads to average over [default: 1]
    --navg <navg>    Average each point over this many frames [default: 50]
    --si             Convert to SI units

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


#def calc_D(xyz, dim, tau, Navg=50):
#    """Diffusion. Input:
#    * xyz: matrix positions (N, dim)
#    * dim: dimension, 1, 2, or 3
#    * tau: time step between frames
#    * Navg: get one D by averaging over this many frames
#    Compute one D by averaging over Navg frames, starting
#    at a specific step, then average over the starting steps."""
#    N = len(xyz)
#    Nl = N - Navg
#    Ds = np.zeros(Nl)
#    if xyz.shape[1] != dim:
#        sys.exit("Matrix of positions does not have correct dimension.")
##    xyz2 = np.sum((xyz - xyz[0])**2, axis=1)
#    for i in range(Nl):
##        Ds[i] = np.sum(xyz2[i:i+Navg]) / (2 * (dim+1) * tau * Navg)
#        Ds[i] = np.sum((xyz[i+Navg] - xyz[i])**2) / (2 * (dim+1) * Navg * tau)
#    return np.average(Ds)


def calc_D_crosses(xyz, dim, tau, Navg=50, cut=1.0, L=10.0):
    """Diffusion. Input:
    * xyz: matrix positions (N, dim)
    * dim: dimension, 1, 2, or 3
    * tau: time step between frames
    * Navg: get one D by averaging over this many frames
    Compute one D by averaging over Navg frames, starting
    at a specific step, then average over the starting steps.
    THINK ABOUT CUTOFF!
    """
    N = len(xyz)
    if xyz.shape[1] != dim:
        sys.exit("Matrix of positions does not have correct dimension.")
    Nl = N - Navg       # Number of successive frames to track down r^2
    Ds = np.zeros(Navg)

    d_xyz = diff_xyz(xyz, L, cut)

    for i in range(Navg):   # Average over steps
        dr2 = sum(np.sum(d_xyz[i:i+Nl], 0)**2)
        Ds[i] = dr2 / (2 * (dim + 1) * Nl * tau)
    return np.average(Ds), np.std(Ds)


def diff_xyz(xyz, L, cut=1.0):
    """Create mat"""
    N, dim = xyz.shape
    d_xyz = np.zeros((N-1, dim))
    for i in range(N-1):
        d_xyz[i] = xyz[i+1] - xyz[i]

#    bins = np.linspace(0, L, 101)
#    Rs = np.sqrt(np.sum(d_xyz**2, 1))
#    plt.hist(Rs, bins=bins)
#    plt.title("Before")
#    plt.show()

    d_xyz[d_xyz > cut] -= L
    d_xyz[d_xyz < -cut] += L
    return d_xyz


def rsq_vs_t(xyz, Navg, L):
    """Plot <r^2> vs time averaged over many particles"""
    N = len(xyz)
    Nl = N - Navg       # walk max this many steps
    R2s = np.zeros(Nl)
    d_xyz = diff_xyz(xyz, L, cut=1.0)
    for l in range(Nl):
        for j in range(Navg):  # averaging
            R2s[l] = sum(np.sum(d_xyz[j:j+l], 0)**2)
        R2s[l] /= Nl
    return R2s


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
    bt = int(args["--bt"])
    check_beadtypes(bt, frames[0])
    Navg = int(args["--navg"])
    Nf = len(frames)
    if Nf < Navg:
        sys.exit("Only %i frames found, at least %i required." % (Nf, Navg))
    dt, freq, L = parse_control()
    tau = dt * freq
    Nb = int(args["--nbeads"])
    cut = 5.0
    if Nb < 1:
        sys.exit("Enter number of beads >= 1.")

    print("===== Diffusivity =====")
    print("Frames: %i | dt: %.3f | steps/frame: %i" % (Nf, dt, freq))
    print("Beads to use: %i | PBC cutoff: %.2f" % (Nb, cut))

    # Reading frames
    xyzs = np.zeros((Nb, 3, Nf))   # (beads, 3, frames)
    for i in range(Nf):
        xyzs[:, :, i] = read_from_top(frames[i], Nb, btype=bt)
    Ds, stdDs = np.zeros(Nb), np.zeros(Nb)

    if args["1d"]:
        axes = {"x": 0, "y": 1, "z": 2}
        ax = args["<ax>"]
        if ax not in axes.keys():
            sys.exit("Enter 'x', 'y', 'z' for axis.")
        ax = axes[ax]
        ax = tuple((ax,))   # matrix slice by tuple returns a np.matrix
        print("1d diffusivity | axis: %i" % ax)
        for i in range(Nb):
            Ds[i], stdDs[i] = calc_D_crosses(xyzs[i, ax, :].T, 1, tau, Navg, cut=cut)

        if args["--plot-rsq"]:
            print("Plotting <r^2> vs t.")
            Nl = Nf - Navg
            t = np.arange(0, Nl * tau, tau)
            R2, R2s = np.zeros(Nl), np.zeros(Nl)
            endRs = np.zeros(Nb)

            ti = time.time()
            for i in range(Nb):
                R2 = rsq_vs_t(xyzs[i, ax, :].T, Navg, L)
                print(R2[-1])
                R2s += R2 / Nb
                endRs[i] = R2[-1]
            R2avg, R2std = np.average(endRs), np.std(endRs)
            tf = time.time()

            print("Average end distance: %.3f | std: %.3f | ratio. %.3f" \
                    % (R2avg, R2std, R2std / R2avg))
            print("Time: %.2f s." % (tf - ti))
            plt.plot(t, R2s)
            plt.show()
            sys.exit()

    if args["2d"]:
        planes = {"xy": (0, 1), "yz": (1, 2), "xz": (0, 2)}
        pl = args["<pl>"]
        if pl not in planes.keys():
            sys.exit("Enter 'xy', 'yz', 'xz' for plane.")
        pl = planes[pl]
        print("2d diffusivity | plane: (%i, %i)" % (pl[0], pl[1]))
        for i in range(Nb):
            Ds[i], stdDs[i] = calc_D_crosses(xyzs[i, pl, :].T, 2, tau, Navg, cut=cut)

        if args["--plot-rsq"]:
            print("Plotting <r^2> vs t.")
            Nl = Nf - Navg
            t = np.arange(0, Nl * tau, tau)
            R2, R2s = np.zeros(Nl), np.zeros(Nl)
            endRs = np.zeros(Nb)

            ti = time.time()
            for i in range(Nb):
                R2 = rsq_vs_t(xyzs[i, pl, :].T, Navg, L)
                print(R2[-1])
                R2s += R2 / Nb
                endRs[i] = R2[-1]
            R2avg, R2std = np.average(endRs), np.std(endRs)
            tf = time.time()

            print("Average end distance: %.3f | std: %.3f | ratio. %.3f" \
                    % (R2avg, R2std, R2std / R2avg))
            print("Time: %.2f s." % (tf - ti))
            plt.plot(t, R2s)
            plt.show()
            sys.exit()

    if args["3d"]:
        if args["--plot-rsq"]:
            print("Plotting <r^2> vs t.")
            Nl = Nf - Navg
            t = np.arange(0, Nl * tau, tau)
            R2, R2s = np.zeros(Nl), np.zeros(Nl)
            endRs = np.zeros(Nb)

            ti = time.time()
            for i in range(Nb):
                R2 = rsq_vs_t(xyzs[i, :, :].T, Navg, L)
                print(R2[-1])
                R2s += R2 / Nb
                endRs[i] = R2[-1]
            R2avg, R2std = np.average(endRs), np.std(endRs)
            tf = time.time()

            print("Average end distance: %.3f | std: %.3f | ratio. %.3f" \
                    % (R2avg, R2std, R2std / R2avg))
            print("Time: %.2f s." % (tf - ti))
            plt.plot(t, R2s)
            plt.show()

        print("Computing 3d (bulk) diffusivity")
        for i in range(Nb):
            Ds[i], stdDs[i] = calc_D_crosses(xyzs[i, :, :].T, 3, tau, Navg, cut=cut)
    
#    print(Ds)
    D = np.average(Ds)
    stdD = np.std(Ds)
    mu_stdD = np.average(stdDs)
    std_stdD = np.std(stdDs)
#    print("Average over: %i frames | D: %.3f | std(D): %.3f | ratio: %.3f" \
#            % (Navg, D, stdD, stdD / D))
    print("Average over: %i frames | D: %.3f | std(D): %.3f" % (Navg, D, stdD))
#    print("mu(std(D)): %.3f | std(std(D): %.3f" % (mu_stdD, std_stdD))

    if args["--si"]:
        print("In SI units (m^2/s):")
        print("Average over: %i frames | D: %.3e | std(D): %.3e" \
                % (Navg, D * r_DPD**2 / tau_DPD, stdD * r_DPD**2 / tau_DPD))
             

