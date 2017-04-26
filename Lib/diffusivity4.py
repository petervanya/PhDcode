#!/usr/bin/env python
"""Usage:
    diffusivity4.py <frames> [--bt <bt> --nb <nb> --si --fit-stop <st>]

Compute mean square distance and diffusivity in linear (Einstein) regime,
i.e. early after the start of measurement.
D = <r^2> / (n t), n = 2 (1d), 4 (2d), 6 (3d)

Options:
    --bt <bt>           Bead type [default: 1]
    --nb <nb>           Number of beads for averaging, or 'all' [default: 1]
    --si                Convert to SI units
    --fit-stop <st>     Fit line from 0 until this time [default: 20]

pv278@cam.ac.uk, 26/04/17
"""
import numpy as np
from numba import jit, float64, int64
import time, os, sys, glob
from docopt import docopt
from dlms_lib import read_xyzfile

AMU = 1.66e-27
m0 = 6 * 18 * AMU
kB = 1.381e-23
T = 300.0
r_DPD = 8.14e-10
tau_DPD = np.sqrt(m0 * r_DPD**2 / (kB * T))


@jit(float64[:, :](float64[:, :], float64[:]), nopython=True)
def diff_xyz(xyz, L):
    """Create matrix of differences. L: vector of size 3"""
    N, dim = xyz.shape
    d_xyz = np.zeros((N-1, dim))
    for i in range(N-1):
        d_xyz[i] = xyz[i+1] - xyz[i]
    d_xyz = d_xyz - (d_xyz > L/2) * L
    d_xyz = d_xyz + (d_xyz < -L/2) * L
    return d_xyz


def gen_rsq(xyz, L):
    """Plot <r^2> vs time averaged over many particles"""
    N = len(xyz)
    rsq = np.zeros((N, 3))
    d_xyz = diff_xyz(xyz, L)
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
    Nball = sum(A[:, 0] == bt)
    if bt not in bts:
        sys.exit("Bead type %i not found in frame." % bt)
    return Nball


def parse_control():
    """Extract timestep from CONTROL file"""
    if not os.path.isfile("CONTROL"):
        sys.exit("CONTROL file not found.")
    f = open("CONTROL", "r").readlines()
    for line in f:
        if "volume" in line:
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
    Nf = len(frames)
    if Nf == 0:
        sys.exit("No frames captured.")
    bt = int(args["--bt"])
    dt, freq, L = parse_control()
    delta_t = dt * freq
    Nball = check_beadtypes(bt, frames[0])
    Nb = args["--nb"]
    if Nb == "all":
        Nb = Nball
    else:
        Nb = int(Nb)
    if Nb > Nball:
        sys.exit("Nb larger than number of beads in frame (%i)." % Nball)

    Rsq = np.zeros((Nf, 3))
    t = np.arange(0, Nf * delta_t, delta_t)
    stop = float(args["--fit-stop"])
    if stop < 0.0 or stop > Nf * delta_t:
        sys.exit("Enter  between 0 and %.1f." % (Nf * delta_t))

    print("===== Diffusivity =====")
    print("Frames: %i | dt: %.3f | steps/frame: %i" % (Nf, dt, freq))
    print("Bead type: %i | Number of beads: %i" % (bt, Nb))
    print("Box size: %s | Time range: %.1f | Fit up to: %.1f" \
            % (L, Nf * delta_t, stop))
    if args["--si"]:
        print("SI units | rc: %.2e | tau: %.2e" % (r_DPD, tau_DPD))

    # Reading frames
    ti = time.time()
    xyzs = np.zeros((Nb, 3, Nf))   # (beads, 3, frames)
    for i in range(Nf):
        xyzs[:, :, i] = read_from_top(frames[i], Nb, btype=bt)
    tf = time.time()
    print("Time reading frames: %.2f s." % (tf - ti))

    # Calculating distances
    ti = time.time()
    for i in range(Nb):
        Rsq += gen_rsq(xyzs[i, :, :].T, L) / Nb
    tf = time.time()
    print("Time computing MSD: %.2f s." % (tf - ti))

    # 1D
    outname = "msd_1d_b%i.out" % bt
    hdr = "t, rsq(x), rsq(y), rsq(z)"
    np.savetxt(outname, np.c_[t, Rsq], fmt="%.6e", header=hdr)

    D = np.zeros(3)
    err = np.zeros(3)
    for i in range(3):
        Cc, Ce = np.polyfit(t[t < stop], Rsq[t < stop, i], 1, cov=True)
        D[i] = Cc[0] / 2.0
        err[i] = np.sqrt(np.diag(Ce))[0] / 2.0
    print("\n1d diffusivity  (x,  y,  z):  %.6f  %.6f  %.6f" % tuple(D))
    print("1d fit error    (x,  y,  z):  %.6f  %.6f  %.6f" % tuple(err))
    print("1d MSD saved in %s." % outname)
    if args["--si"]:
        print("1d in SI:  %.6e  %.6e  %.6e" % tuple(D * r_DPD**2 / tau_DPD))

    # 2D
    Rsq_2d = np.zeros((Nf, 3))
    Rsq_2d[:, 0] = Rsq[:, 0] + Rsq[:, 1]  # xy
    Rsq_2d[:, 1] = Rsq[:, 1] + Rsq[:, 2]  # yz
    Rsq_2d[:, 2] = Rsq[:, 0] + Rsq[:, 2]  # xz

    outname = "msd_2d_b%i.out" % bt
    hdr = "t, rsq(xy), rsq(yz), rsq(xz)"
    np.savetxt(outname, np.c_[t, Rsq_2d], fmt="%.6e", header=hdr)

    D, err = np.zeros(3), np.zeros(3)
    for i in range(3):
        Cc, Ce = np.polyfit(t[t < stop], Rsq_2d[t < stop, i], 1, cov=True)
        D[i] = Cc[0] / 4.0
        err[i] = np.sqrt(np.diag(Ce))[0] / 4.0
    print("\n2d diffusivity (xy, yz, xz):  %.6f  %.6f  %.6f" % tuple(D))
    print("2d fit error   (xy, yz, xz):  %.6f  %.6f  %.6f" % tuple(err))
    print("2d MSD saved in %s." % outname)
    if args["--si"]:
        print("2d in SI:  %.6e  %.6e  %.6e" % tuple(D * r_DPD**2 / tau_DPD))

    # 3D
    Rsq_3d = np.sum(Rsq, 1)
    outname = "msd_3d_b%i.out" % bt
    np.savetxt(outname, np.c_[t, Rsq_3d], fmt="%.6e")

    Cc, Ce = np.polyfit(t[t < stop], Rsq_3d[t < stop], 1, cov=True)
    D = Cc[0] / 6.0
    err = np.sqrt(np.diag(Ce))[0] / 6.0
    print("\n3d diffusivity: %.6f" % D)
    print("3d fit error  : %.6f" % err)
    print("3d MSD saved in %s." % outname)
    if args["--si"]:
        print("3d in SI:  %.6e" % (D * r_DPD**2 / tau_DPD))


