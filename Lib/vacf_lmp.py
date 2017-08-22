#!/usr/bin/env python
"""Usage:
    vacf_lmp.py <frames> [--bt <bt> --nb <nb> --si --start <st>]
                         [--freq <f> --dt <dt>]

Compute velocity autocorrelation function in LAMMPS.
Cv(t) = <v(t) . v(0)>

Options:
    --bt <bt>      Bead type [default: 1]
    --nb <nb>      Number of beads for averaging, or 'all' [default: all]
    --si           Convert to SI units
    --start <st>   Fraction of frames after which to start fitting a line
                   over rsq vs t [default: 0.666]
    --dt <dt>      Timestep [default: 0.05]
    --freq <f>     Dump frequency [default: 1]

pv278@cam.ac.uk, 21/08/17
"""
import numpy as np
from scipy.integrate import simps
import sys, glob, time
from docopt import docopt


def read_velfile(outfile):
    try:
        A = np.loadtxt(outfile, skiprows=9)
    except FileNotFoundError:
        sys.exit("File %s not found." % outfile)
    return A


def check_beadtypes(bt, frame):
    """Check if beadtype present in the xyz frame"""
    A = read_velfile(frame)
    bts = set(A[:, 0].astype(int))
    Nbf = sum(A[:, 0] == bt)
    if bt not in bts:
        sys.exit("Bead type %i not found in frame." % bt)
    return Nbf


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


if __name__ == "__main__":
    args = docopt(__doc__)
    frames = sorted(glob.glob(args["<frames>"]))
    if len(frames) == 0:
        sys.exit("0 frames captured.")
    Nf = len(frames)
    bt = int(args["--bt"])
    dt = float(args["--dt"])
    freq = int(args["--freq"])
    tau = dt * freq
    Nbf = check_beadtypes(bt, frames[0])
    Nb = args["--nb"]
    if Nb == "all":
        Nb = Nbf
    else:
        Nb = int(Nb)
    if Nb > Nbf:
        sys.exit("Nb larger than number of beads in frame (%i)." % Nbf)

    print("===== Velocity autocorrelation function =====")
    print("Frames: %i | dt: %.3f | steps/frame: %i" % (Nf, dt, freq))
    print("Bead type: %i | Beads: %i" % (bt, Nb))

    xyzs = np.zeros((Nb, 3, Nf))   # (beads, dim, frames)
    for i in range(Nf):
        xyzs[:, :, i] = read_from_top(frames[i], Nb, btype=bt)
    Cvv = np.zeros((Nf, 3))
    t = np.arange(0, Nf, 1) * tau

    ti = time.time()
    for i in range(Nf):
        Cvv[i] = np.sum(xyzs[:, :, 0] * xyzs[:, :, i], 0) / Nb
    tf = time.time()
    print("Time: %.2f s." % (tf - ti))

    outname = "cvv_b%i.out" % bt
    hdr = "t, cvv(x), cvv(y), cvv(z)"
    np.savetxt(outname, np.c_[t, Cvv], fmt="%.6e", header=hdr)

    D = np.array([simps(Cvv[:, i], t) for i in range(3)])

    print("1d diffusivity  (x,  y,  z):  %.6f  %.6f  %.6f" % tuple(D))
    print("3d diffusivity: %.6f" % (np.sum(D) / 3.0))


