#!/usr/bin/env python
"""Usage:
    dlms_binmixt.py [--f <f> --L <L> --rho <rho> --xyz --vel <T>]
                    [--aii <aii> --da <da>]

Generate a binary mixture (A/B particles) CONFIG and FIELD file 
for DL_MESO package.

Options:
    --f <f>        Fraction of A partcicles [default: 0.5]
    --L <L>        Box size [default: 10.0]
    --rho <rho>    Density of particles [default: 3.0]
    --aii <aii>    Same bead interaction [default: 25]
    --da <da>      Excess interaction between A and B [default: 5]
    --xyz          Create xyz file
    --vel <T>      Initialise velocities at given temperature T

pv278@cam.ac.uk, 02/06/16, modified 04/10/16
"""
import numpy as np
from lmp_lib import save_xyzfile
from dlms_lib import save_config, inter2str, species2str2
from docopt import docopt
import sys


if __name__ == "__main__":
    args = docopt(__doc__)
    np.random.seed(1234)
    f = float(args["--f"])
    L = float(args["--L"])
    rho = float(args["--rho"])
    aii = float(args["--aii"])
    da = float(args["--da"])
    N = int(rho * L**3) // 10 * 10     # round to 10
    if f < 0.0 or f > 1.0:
        sys.exit("Fraction f of A beads must be between 0 and 1.")
    print("N: %i | L: %.1f | rho: %.1f | fraction: %.2f" % (N, L, rho, f))

    NA = int(f * N)
    NB = N - NA
    names = ["A"] * NA + ["B"] * NB
    xyz = np.random.rand(N, 3) * L

    levcfg = 0
    if args["--vel"]:
        T = float(args["--vel"])
        print("Initialising velocities, temperature: %.1f" % T)
        vel = np.random.randn(N, 3)*T
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
#    save_config("CONFIG", names, xyz)   # use function form dlms_lib.py

    # ===== generate FIELD file with bead species and interactions
    rc, gamma = 1.0, 4.5
    a_ij = {}
    a_ij["A A"] = [aii, rc, gamma]
    a_ij["B B"] = [aii, rc, gamma]
    a_ij["A B"] = [aii + da, rc, gamma]
    field_string = "bla\n\n" +\
                   species2str2({"A": NA, "B": NB}) +\
                   inter2str(a_ij) +\
                   "close\n"
    open("FIELD", "w").write(field_string)
    print("FIELD file saved.")


    if args["--xyz"]:
        fname = "CONFIG_INIT.xyz"
        names = [1] * NA + [2] * NB
        save_xyzfile(fname, np.hstack((np.matrix(names).T, xyz)) )
        print("xyz file written in", fname)


