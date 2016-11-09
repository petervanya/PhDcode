#!/usr/bin/env python
"""Usage: integrate_profile.py <profile> <d> [--L <L>]

[AD HOC] Load 1d water profile and integrate 
volume of water in polymer and in electrodes.

Arguments:
    <file>    Water profile, columns [r, f]
    <d>       Slab width in nm

Options:
    --L <L>   Box size in DPD units [default: 40]

09/11/16
"""
import numpy as np
from scipy.integrate import simps
import sys
from docopt import docopt

rc = 8.14e-10

if __name__ == "__main__":
    args = docopt(__doc__)
    L = float(args["--L"])
    d_nm = float(args["<d>"])
    d = d_nm * 1e-9 / rc
    A = np.loadtxt(args["<profile>"])
    r, f = A[:, 0], A[:, 1]
    if d < 0.0 or d > L:
        sys.exit("Slab width larger than box size.")
    print("===== Integrating water profile =====")
    print("L: %.2f | slab width:  %.2f (%.2f nm)" % (L, d, d_nm))

    dr = r[1] - r[0]
    re1 = r[r < (L-d)/2]
    re2 = r[r > (L+d)/2]
    rm = r[(r >= (L-d)/2) & (r <= (L+d)/2)]

    fe1 = f[r < (L-d)/2]
    fe2 = f[r > (L+d)/2]
    fm = f[(r >= (L-d)/2) & (r <= (L+d)/2)]

    water_mat = simps(fm, dx=dr)
    water_elec = simps(fe1, dx=dr) + simps(fe2, dx=dr)
    water_tot = simps(f, dx=dr)

    print("Total water: %.2f" % water_tot)
    print("Electrodes: %.2f | Matrix: %.2f | mat / el: %.2f" % \
            (water_elec, water_mat, water_mat / water_elec))
    
#    water_mat = np.sum(fm) * dr
#    water_elec = (np.sum(fe1) + np.sum(fe2)) * dr
#    water_tot = np.sum(f) * dr
#
#    print("Naive quadrature | Total water: %.2f" % water_tot)
#    print("Electrodes: %.2f | Matrix: %.2f | mat / el: %.2f" % \
#            (water_elec, water_mat, water_mat / water_elec))


