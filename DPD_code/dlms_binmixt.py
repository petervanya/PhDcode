#!/usr/bin/env python
"""Usage:
    dlms_binmixt.py [--f <f> --L <L> --rho <rho> --xyz --vel <T>]
                    [--aii <aii> --da <da> --loc --gamma <g>]

Generate a binary mixture (A/B particles) CONFIG and FIELD file 
for DL_MESO package.

Options:
    --f <f>        Fraction of A partcicles [default: 0.5]
    --L <L>        Box size, 'L' or 'Lx Ly Lz' [default: 10.0]
    --rho <rho>    Density of particles [default: 3.0]
    --aii <aii>    Same bead interaction [default: 25]
    --da <da>      Excess interaction between A and B [default: 5]
    --xyz          Create xyz file
    --vel <T>      Initialise velocities at given temperature T
    --loc          Localise A beads in left and B beads in right side of box
    --gamma <g>    Friction [default: 4.5]

pv278@cam.ac.uk, 02/06/16, modified 15/02/17
"""
import numpy as np
from dlms_lib import inter2str, species2str, save_xyzfile
import sys
from docopt import docopt


if __name__ == "__main__":
    args = docopt(__doc__)
    np.random.seed(1234)
    f = float(args["--f"])
    rho = float(args["--rho"])
    aii = float(args["--aii"])
    da = float(args["--da"])
    s = args["--L"].split()
    if len(s) == 1:
        L = float(s[0]) * np.ones(3)
    elif len(s) == 3:
        L = np.array(s).astype(float)
    else:
        sys.exit("<L> should have size 1 or 3.")
    N = int(rho * np.prod(L))
    if f < 0.0 or f > 1.0:
        sys.exit("Fraction f of A beads must be between 0 and 1.")

    print("===== DPD binary mixture for DL_MESO =====")
    print("N: %i | rho: %.1f | fraction: %.2f" % (N, rho, f))
    print("Box: %s" % str(L))

    NA = int(f * N)
    NB = N - NA
    names = ["A"] * NA + ["B"] * NB
    xyz = np.random.rand(N, 3) * L
    if args["--loc"]:
        xyz[:NA, 0] = np.random.rand(NA) * f * L[0]
        xyz[NA:, 0] = np.random.rand(NB) * (1 - f) * L[0] + f * L[0]

    levcfg = 0
    if args["--vel"]:
        T = float(args["--vel"])
        print("Initialising velocities, temperature: %.1f" % T)
        vel = np.random.randn(N, 3) * T
        vel -= np.sum(vel, 0) / N
        levcfg = 1

    conf_str = "bla\n" + "0\t" + str(levcfg) + "\n"
    for i in range(N):
        conf_str += "%s        %i\n" % (names[i], (i+1))
        conf_str += "    %.10f    %.10f    %.10f\n" % \
                (xyz[i, 0], xyz[i, 1], xyz[i, 2])
        if levcfg == 1:
            conf_str += "    %.10f    %.10f    %.10f\n" % \
                    (vel[i, 0], vel[i, 1], vel[i, 2])

    fname = "CONFIG"
    open(fname, "w").write(conf_str)
    print("Initial configuration written in %s" % fname)

    rc, gamma = 1.0, float(args["--gamma"])
    a_ij = {}
    a_ij["A A"] = [aii, rc, gamma]
    a_ij["B B"] = [aii, rc, gamma]
    a_ij["A B"] = [aii + da, rc, gamma]
    field_string = "bla\n\n" + \
            species2str({"A": NA, "B": NB}) + \
            inter2str(a_ij) + \
            "close\n"
    open("FIELD", "w").write(field_string)
    print("FIELD file saved.")

    if args["--xyz"]:
        fname = "CONFIG_INIT.xyz"
        names = [1] * NA + [2] * NB
        save_xyzfile(fname, np.hstack((np.matrix(names).T, xyz)))
        print("xyz file written in", fname)


