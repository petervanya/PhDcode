#!/usr/bin/env python
"""Usage:
    dlms_binmixt_mdpd.py [--f <f> --L <L> --rho <rho> --xyz --vel <T>]
                         [--A <A> --dA <dA> --B <B> --dB <dB>]

Generate a binary mixture (A/B particles) CONFIG and FIELD file 
for DL_MESO package.

Options:
    --L <L>        Box size, either 'L' or 'Lx Ly Lz' [default: 10.0]
    --f <f>        Fraction of A partcicles [default: 0.5]
    --rho <rho>    Density of particles [default: 3.0]
    --A <A>        Same bead interaction [default: -25]
    --dA <dA>      Excess interaction between A and B [default: 5]
    --B <B>        Same bead interaction [default: 25]
    --dB <dB>      Excess interaction between A and B [default: 0]
    --xyz          Create xyz file
    --vel <T>      Initialise velocities at given temperature T

pv278@cam.ac.uk, 03/02/17
"""
import numpy as np
from dlms_lib import save_config, save_xyzfile, inter2str_mdpd, species2str
import sys
from docopt import docopt


if __name__ == "__main__":
    args = docopt(__doc__)
    np.random.seed(1234)
    f = float(args["--f"])
    s = args["--L"].split()
    if len(s) == 1:
        L = float(s[0]) * np.ones(3)
    elif len(s) == 3:
        L = np.array(s).astype(float)
    else:
        sys.exit("<L> should have size 1 or 3.")
    rho = float(args["--rho"])
    A = float(args["--A"])
    dA = float(args["--dA"])
    B = float(args["--B"])
    dB = float(args["--dB"])
    N = int(rho * np.prod(L)) // 10 * 10     # round to 10
    if f < 0.0 or f > 1.0:
        sys.exit("Fraction f of A beads must be between 0 and 1.")

    print("===== Many-body DPD binary mixture =====")
    print("N: %i | rho: %.1f | fraction: %.2f" % (N, rho, f))
    print("Box:", L)

    NA = int(f * N)
    NB = N - NA
    names = ["A"] * NA + ["B"] * NB
    xyz = np.random.rand(N, 3) * L
    save_config("CONFIG", names, xyz, np.diag(L))

    # ===== generate FIELD file with bead species and interactions
    rc, gamma = 1.0, 4.5
    a_ij = {}
    a_ij["A A"] = [A, B, 0, 0, 0, rc, gamma]
    a_ij["B B"] = [A, B, 0, 0, 0, rc, gamma]
    a_ij["A B"] = [A + dA, B + dB, 0, 0, 0, rc, gamma]
    field_string = "bla\n\n" +\
                   species2str({"A": NA, "B": NB}) +\
                   inter2str_mdpd(a_ij) +\
                   "close\n"
    open("FIELD", "w").write(field_string)
    print("FIELD file saved.")


    if args["--xyz"]:
        fname = "CONFIG_INIT.xyz"
        names = [1] * NA + [2] * NB
        save_xyzfile(fname, np.hstack((np.matrix(names).T, xyz)) )
        print("xyz file written in", fname)


