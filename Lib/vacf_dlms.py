#!/usr/bin/env python
"""Usage:
    vacf.py <frames> [--bt <bt> --nb <nb> --si --mb <mb> --start <st>]

Compute velocity autocorrelation function in DL_MESO.
Cv(t) = <v(t) . v(0)>


Options:
    --bt <bt>      Bead type [default: 1]
    --nb <nb>      Number of beads for averaging, or 'all' [default: all]
    --si           Convert to SI units
    --mb <mb>      Rescale by number of molecules in a bead [default: 1]
    --start <st>   Fraction of frames after which to start fitting a line
                   over rsq vs t [default: 0.666]

pv278@cam.ac.uk, 12/04/17, modified 22/08/17
"""
import numpy as np
from scipy.integrate import simps
import sys, glob, time
from docopt import docopt
from dlms_lib import read_xyzfile


AMU = 1.66e-27
m0 = 6 * 18 * AMU
kB = 1.381e-23
T = 300.0
r_DPD = 8.14e-10
tau_DPD = np.sqrt(m0 * r_DPD**2 / (kB * T))


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

    xyzs = np.zeros((Nb, 3, Nf))   # (beads, dim, frames)
    for i in range(Nf):
        xyzs[:, :, i] = read_from_top(frames[i], Nb, btype=bt)
    Cvv = np.zeros((Nf, 3))
    t = np.arange(0, Nf, 1) * tau

    ti = time.time()
    for i in range(Nb):
        Cvv[i] = np.sum(xyzs[:, :, 0] * xyzs[:, :, i], 0) / Nb
    tf = time.time()
    print("Time: %.2f s." % (tf - ti))

    outname = "cvv_b%i.out" % bt
    hdr = "t, cvv(x), cvv(y), cvv(z)"
    np.savetxt(outname, np.c_[t, Cvv], fmt="%.6e", header=hdr)

    D = np.array([simps(Cvv[:, i], t) for i in range(3)])

    print("1d diffusivity  (x,  y,  z):  %.6f  %.6f  %.6f" % tuple(D))
    print("3d diffusivity: %.6f" % (np.sum(D) / 3.0))


