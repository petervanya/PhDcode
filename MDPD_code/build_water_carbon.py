#!/usr/bin/env python
"""Usage:
    build_water_carbon.py [--L <L> --hc <hc> --lw <lw>]
                          [--A <A> --B <B> --rho <rho> --chi <chi>]
                          [--gamma <g> --xyz --grav]

Build a water cube on top of a carbon slab of width wc.
+--------+
|  www   |
|  www   |
|CCCCCCCC|
+--------+
Water-carbon chi-param taken from Wu, EES, 2008.

Options:
    --L <L>       Box size [default: 10]
    --hc <hc>     Height of carbon [default: 3]
    --lw <lw>     Side of a water cube [default: 3]
    --rho <rho>   Density [default: 6]
    --A <A>       Attraction [default: -80.4]
    --B <B>       Repulsion [default: 51.7]
    --chi <chi>   Carbon-water chi-param [default: 3.7653]
    --gamma <g>   Friction [default: 15]
    --grav        Include gravity

pv278@cam.ac.uk, 07/12/17
"""
import numpy as np
from dlms_lib import save_config, save_xyzfile, inter2str_mdpd, \
        species2str2, chi2da_mdpd
import sys
from docopt import docopt


if __name__ == "__main__":
    args = docopt(__doc__)
    np.random.seed(1234)
    rho = float(args["--rho"])
    A = float(args["--A"])
    B = float(args["--B"])
    chi = float(args["--chi"])
    hc = float(args["--hc"])
    lw = float(args["--lw"])
    s = args["--L"].split()
    if len(s) == 1:
        L = float(s[0]) * np.ones(3)
    elif len(s) == 3:
        L = np.array(s).astype(float)
    else:
        sys.exit("<L> should have size 1 or 3.")

    Nc = int(rho * 2 * L[0] * L[1] * hc)
    Nw = int(rho * lw**3)

    print("===== Many-body DPD water cube on carbon slab =====")
    print("rho: %.1f | slab width: %.3f | water size: %.3f" \
            % (rho, hc, lw))
    print("C beads: %i | W beads: %i | Box: %s" % (Nc, Nw, L))

    xyzc = np.random.rand(Nc, 3) * L
    xyzc[:, 2] *= hc / L[2]
    xyzw = np.random.rand(Nw, 3) * lw
    xyzw += [(L[0]-lw)/2, (L[1]-lw)/2, hc]
    names = ["C"] * Nc + ["W"] * Nw
    xyz = np.r_[xyzc, xyzw]

    conf_str = "bla\n0    0\n"
    for i in range(Nc+Nw):
        conf_str += "%s        %i\n" % (names[i], (i+1))
        conf_str += "    %.10f    %.10f    %.10f\n" \
            % (xyz[i, 0], xyz[i, 1], xyz[i, 2])

    fname = "CONFIG"
    open(fname, "w").write(conf_str)
    print("Initial configuration written in %s" % fname)

    # ===== FIELD file with bead species and interactions
    rc = 1.0
    gamma = float(args["--gamma"])
    dA = chi2da_mdpd(A, B) * chi

    if args["--grav"]:
        AMU = 1.66e-27
        m = 6 * 18 * AMU
        kT = 1.381e-23 * 300.0
        g = -10.0 * m / kT

    a_ij = {}
    a_ij["C C"] = [A, B, 0, 0, 0, rc, gamma]
    a_ij["W W"] = [A, B, 0, 0, 0, rc, gamma]
    a_ij["C W"] = [A + dA, B, 0, 0, 0, rc, gamma]
    field_string = "bla\n\n" + \
            species2str2({"C": [Nc, 1], "W": [Nw, 0]}) + \
            inter2str_mdpd(a_ij)
    if args["--grav"]:
        field_string += "external\ngrav 0.0 0.0 %.6f\n\n" % g
    field_string += "close\n"
    open("FIELD", "w").write(field_string)
    print("FIELD file saved.")

    if args["--xyz"]:
        fname = "CONFIG_INIT.xyz"
        names = [1] * Nc + [2] * Nw
        save_xyzfile(fname, np.c_[names, xyz])
        print("xyz file written in %s." % fname)



