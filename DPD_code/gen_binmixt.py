#!/usr/bin/env python
"""Usage:
    gen_binmixt.py [--f <f> --L <L> --rho <rho> --xyz --vel <T> --loc]

Create a LAMMPS data file for A/B binary mixture with fraction f of A beads.
LAMMPS atom_style: atomic.

Options:
    --f <f>        Fraction of A beads [default: 0.5]
    --L <L>        Box size [default: 10]
    --rho <rho>    Density of the system [default: 3.0]
    --xyz          Create xyz file
    --vel <T>      Initialise velocities at given temperature T
    --loc          Localise A beads in left and B beads in right side of box

pv278@cam.ac.uk, 13/06/16
"""
import numpy as np
from lmp_lib import save_xyzfile
from docopt import docopt
import sys


def header2str(N, L):
    """Generate LAMMPS header"""
    s = "#binmixt\n"
    s += str(N) + " atoms\n"
    s += "2 atom types\n"
    s += "\n"
    s += "0.0 " + str(L[0]) + " xlo xhi\n"
    s += "0.0 " + str(L[1]) + " ylo yhi\n"
    s += "0.0 " + str(L[2]) + " zlo zhi\n\n"
    return s


def mass2str(m1, m2):
    s = "Masses\n\n"
    s += "1 " + str(m1) + "\n"
    s += "2 " + str(m2) + "\n\n"
    return s


def atoms2str(names, xyz):
    """Convert atomic matrix to string"""
    N = len(xyz)
    s = "Atoms\n\n"
    for i in range(N):
        s += "%i\t%i\t%f\t%f\t%f\n" % \
             (i+1, names[i], xyz[i, 0], xyz[i, 1], xyz[i, 2])
    s += "\n"
    return s


def vel2str(vel):
    """Convert atomic matrix to string"""
    N = len(vel)
    s = "Velocities\n\n"
    for i in range(N):
        s += "%i\t%f\t%f\t%f\n" % (i+1, vel[i, 0], vel[i, 1], vel[i, 2])
    s += "\n"
    return s


if __name__ == "__main__":
    args = docopt(__doc__)
    np.random.seed(1234)
    f = float(args["--f"])
    rho = float(args["--rho"])
    s = args["--L"].split()
    if len(s) == 1:
        L = float(s[0]) * np.ones(3)
    elif len(s) == 3:
        L = np.array(s).astype(float)
    else:
        sys.exit("<L> should have size 1 or 3.")
    N = int(rho * np.prod(L)) // 10 * 10     # round to 10
    if f < 0.0 or f > 1.0:
        print("Fraction f of A beads must be between 0 and 1.")
        sys.exit()

    print("=== LAMMPS data file for binary mixture ====")
    print("N: %i | rho: %.1f | fraction: %.2f" % (N, rho, f))
    print("Box: %s" % str(L))

    NA = int(f * N)
    NB = N - NA
    names = [1] * NA + [2] * NB
    xyz = np.random.rand(N, 3) * L
    if args["--loc"]:
        xyz[:NA, 0] = np.random.rand(NA) * f * L[0]
        xyz[NA:, 0] = np.random.rand(NB) * (1 - f) * L[0] + f * L[0]

    header = header2str(N, L)
    final_string = header + \
                   mass2str(1.0, 1.0) + \
                   atoms2str(names, xyz)

    if args["--vel"]:
        T = float(args["--vel"])
        print("Initialising velocities, temperature: %.1f" % T)
        vel = np.random.randn(N, 3) * T
        vel -= np.sum(vel, 0) / N
        final_string += vel2str(vel)

    fname = "binmixt.data"
    open(fname, "w").write(final_string)
    print("Data file written in", fname)

    if args["--xyz"]:
        fname = "binmixt.xyz"
        save_xyzfile(fname, np.c_[names, xyz])
        print("xyz file written in", fname)


