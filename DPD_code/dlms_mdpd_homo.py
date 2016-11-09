#!/usr/bin/env python
"""Usage:
    dlms_mdpd_homo.py [--L <L> --rho <rho> --xyz]
                      [--A <A> --B <B> --rd <rd>]
                      [--dt <dt> --steps <ns> --traj-after <ta>]

Generate a many-body DPD homogenous mixture.
CONFIG, FIELD and CONTROL files.
for DL_MESO package.

Options:
    --L <L>             Box size [default: 10.0]
    --rho <rho>         Density of particles [default: 3.0]
    --A <A>             Interaction parameter A [default: 25]
    --B <B>             Interaction parameter B [default: 25]
    --rd <rd>           Second cutoff, less than 1 [default: 0.75]
    --dt <dt>           Time step [default: 0.03]
    --steps <ns>        Total number of steps [default: 80000]
    --traj-after <ta>   Start taking frames after n steps [default: 70000]
    --xyz               Create xyz file

pv278@cam.ac.uk, 28/10/16
"""
import numpy as np
from lmp_lib import save_xyzfile
from dlms_lib import save_config, inter2str_mdpd, species2str2, gen_control
from docopt import docopt
import sys


if __name__ == "__main__":
    args = docopt(__doc__)
    np.random.seed(1234)
    L = float(args["--L"])
    rho = float(args["--rho"])
    A = float(args["--A"])
    B = float(args["--B"])
    rd = float(args["--rd"])
    if rd < 0.0 or rd > 1.0:
        sys.exit("Enter rd between 0 and 1.")
    N = int(rho * L**3) // 10 * 10     # round to 10
    print("===== MDPD homogenous mixture =====")
    print("N: %i | L: %.1f | rho: %.1f" % (N, L, rho))
    print("Interaction. A: %.2f | B: %.2f | rd: %.2f" % (A, B, rd))

    names = ["A"] * N
    xyz = np.random.rand(N, 3) * L

    conf_str = "bla\n0    0\n"
    for i in range(N):
        conf_str += "%s        %i\n" % (names[i], (i+1))
        conf_str += "    %.10f    %.10f    %.10f\n" \
            % (xyz[i, 0], xyz[i, 1], xyz[i, 2])

    fname = "CONFIG"
    open(fname, "w").write(conf_str)
    print("Initial configuration written in %s" % fname)

    # ===== FIELD file with bead species and interactions
    rc, gamma = 1.0, 4.5
    a_ij = {}
    a_ij["A A"] = [A, B, 0.0, 0.0, 0.0, rc, gamma]
    field_string = "bla\n\n" +\
                   species2str2({"A": N}) +\
                   inter2str_mdpd(a_ij) +\
                   "close\n"
    open("FIELD", "w").write(field_string)
    print("FIELD file saved.")

    # ==== CONTROL file
    steps = int(args["--steps"])
    ta = int(args["--traj-after"])
    if ta > steps:
        sys.exit("<ta> must be less than number of steps.")
    dt = float(args["--dt"])
    gen_control(L, dt, steps, traj_after=ta, method="mdpd", rd=rd)

    if args["--xyz"]:
        fname = "CONFIG_INIT.xyz"
        names = [1] * N
        save_xyzfile(fname, np.hstack((np.matrix(names).T, xyz)) )
        print("xyz file written in %s." % fname)


