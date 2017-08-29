#!/usr/bin/env python
"""Usage:
    vacf_dlms.py <frames> [--bt <bt> --nb <nb> --si --mb <mb>]

Compute velocity autocorrelation function in DL_MESO.
Cv(t) = <v(t) . v(0)>


Options:
    --bt <bt>      Bead type [default: 1]
    --nb <nb>      Number of beads for averaging, or 'all' [default: all]
    --si           Convert to SI units
    --mb <mb>      Rescale by number of molecules in a bead [default: 1]

pv278@cam.ac.uk, 12/04/17, modified 22/08/17
"""
import numpy as np
from scipy.integrate import simps
import os, sys, glob, time
from docopt import docopt
from dlms_lib import read_xyzfile2


def check_beadtypes(bt, frame):
    """Check if beadtype present in the xyz frame"""
    nm, xyz = read_xyzfile2(frame)
    if bt not in set(nm):
        sys.exit("Bead type %i not found in frame." % bt)
    Nbf = sum(nm == bt)
    return Nbf


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


#def gen_cvv(xyz):
#    Nt = len(xyz)     # Number of timesteps
#    Cvv = np.zeros((Nt-1, 3))
#    for i in range(Nt-1):
#        Cvv[i] = xyz[0] * xyz[i+1]
#    return Cvv


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

    print("===== Velocity autocorrelation function =====")
    print("Frames: %i | dt: %.3f | steps/frame: %i" % (Nf, dt, freq))
    print("Bead type: %i | Beads: %i | Molecules/bead: %i" % (bt, Nb, Mb))

    xyzs = np.zeros((Nbf, 3, Nf))   # (beads, dim, frames)
    ti = time.time()
    for i in range(Nf):
        xyzs[:, :, i] = read_from_top(frames[i], Nb, btype=bt)
    tf = time.time()
    print("Time reading frames: %.2f s." % (tf - ti))
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
    print("VACF saved in %s." % outname)

    D = np.array([simps(Cvv[:, i], t) for i in range(3)])
    D_2d = np.array([D[0] + D[1], D[1] + D[2], D[0] + D[2]]) / 2.0
    D_3d = np.sum(D) / 3.0

    print("1d diffusivity  (x,  y,  z):  %.6f  %.6f  %.6f" % tuple(D))
    print("2d diffusivity (xy, yz, xz):  %.6f  %.6f  %.6f" % tuple(D_2d))
    print("3d diffusivity: %.6f" % D_3d)


