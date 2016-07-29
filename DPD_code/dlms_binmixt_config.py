#!/usr/bin/env python
"""Usage:
    dlms_binmixt_config.py [--f <f> --L <L> --rho <rho> --xyz --vel <T>]

Generate binary mixture (A/B particles) initial configuration file
for DL_MESO package.

Options:
    --f <f>        Fraction of A partcicles [default: 0.5]
    --L <L>        Box size [default: 10.0]
    --rho <rho>    Density of particles [default: 3.0]
    --xyz          Create xyz file
    --vel <T>      Initialise velocities at given temperature T

pv278@cam.ac.uk, 02/06/16
"""
import numpy as np
import lmp_lib as ll
from docopt import docopt
import sys


if __name__ == "__main__":
    args = docopt(__doc__)
    np.random.seed(1234)
    f = float(args["--f"])
    L = float(args["--L"])
    rho = float(args["--rho"])
    N = int(rho * L**3)
    if f < 0.0 or f > 1.0:
        print("Fraction f of A beads must be between 0 and 1.")
        sys.exit()
    print("L: %.1f | Density: %.1f | A beads fraction: %.1f" % (L, rho, f))

    NA, NB = int(f*N), int((1-f)*N)
    names = ["A"]*NA + ["B"]*NB
    xyz = np.random.rand(N, 3)*L
    levcfg = 0
    if args["--vel"]:
        T = float(args["--vel"])
        print("Initialising velocities, temperature: %.1f" % T)
        vel = np.random.randn(N, 3)*T
        levcfg = 1

    conf_str = "pokus\n" + "0\t" + str(levcfg) + "\n"
    for i in range(N):
        conf_str += "%s        %i\n" % (names[i], (i+1))
        conf_str += "    %.10f    %.10f    %.10f\n" % (xyz[i, 0], xyz[i, 1], xyz[i, 2])
        if levcfg == 1:
            conf_str += "    %.10f    %.10f    %.10f\n" % (vel[i, 0], vel[i, 1], vel[i, 2])

    fname = "CONFIG"
    open(fname, "w").write(conf_str)
    print("Initial configuration written in %s" % fname)

    if args["--xyz"]:
        fname = "CONFIG_INIT.xyz"
        names = [1]*NA + [2]*NB
        ll.save_xyzfile(fname, np.hstack((np.matrix(names).T, xyz)) )
        print("xyz file written in", fname)


