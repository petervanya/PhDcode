#!/usr/bin/env python
"""Usage:
    predict_dpd_coeff.py <solvent> <peaksfile> [--rho <rho>]

Read 'dpd_coeffs_rho<n>.out' file (fit: a*(x-x0)**b)
for a specific density, and 'max.out' file with MD peaks,
and predict the interaction coeffs and bead sizes.

Sources
=======
Hanbook of chemistry and physics
http://pubs.acs.org/doi/pdf/10.1021/je060415l
http://www.engineeringtoolbox.com/bulk-modulus-elasticity-d_585.html
http://aip.scitation.org/doi/pdf/10.1063/1.1679903

Options:
    --rho <rho>    DPD density [default: 3]

08/04/17
"""
import numpy as np
import pandas as pd
import sys
from docopt import docopt


dpd_peak_pos = {}
dpd_peak_pos[3] = [9.28038e-01, 6.74471e-01, -2.38334e-02]

dpd_peak_max = {}
dpd_peak_max[3] = [6.86478e-01, -1.21659e+01, 1.44429e-01]


T = 300.0
kB = 1.38e-23
NA = 6.022e23

subst = pd.DataFrame(columns = ["rho", "Mm", "sol", "kappaT"])
subst.loc["SPC"] = [997.0, 18.02e-3, 47.8, 4.504e-10]
subst.loc["TIP3P_EW"] = [997.0, 18.02e-3, 47.8, 4.504e-10]
subst.loc["TIP4P_2005"] = [997.0, 18.02e-3, 47.8, 4.504e-10]
subst.loc["Methanol"] = [792.0, 32.04e-3, 29.7, 1.273e-9]
subst.loc["Ethanol"] = [789.0, 46.07e-3, 26.5, 1.185e-9]
subst.loc["Acetone"] = [791.0, 58.08e-3, 19.9, 1.262e-9]   # at 20 C
subst.loc["SPC_T320"] = [997.0, 18.02e-3, 47.8, 4.170e-10]
subst.loc["SPC_T340"] = [997.0, 18.02e-3, 47.8, 4.500e-10]
subst.loc["SPC_T360"] = [997.0, 18.02e-3, 47.8, 4.700e-10]

subst["n"] = subst.rho * NA / subst.Mm
subst["inv_kappa"] =  1.0 / (subst.n * kB * T * subst.kappaT)
subst["a"] = (subst.inv_kappa - 1.0) / (2 * 0.101)


def fit_func(x, a, x0, b):
    return a * (x - x0)**b


def fit_func_inv(y, a, x0, b):
    return (y / a)**(1.0/b) + x0


args = docopt(__doc__)
sol = args["<solvent>"].rstrip("/")
rho = int(args["--rho"])
if rho not in dpd_peak_max.keys():
    sys.exit("Fit for density %i is not present." % rho)

try:
    peaks_file = args["<peaksfile>"]
    B = np.loadtxt(peaks_file)
    Nm, pos, peaks = B[:, 0], B[:, 1], B[:, 2]
except FileNotFoundError:
    sys.exit("File %s not found." % peaks_file)

try:
    pred_as = fit_func_inv(peaks, *dpd_peak_max[rho])
    pred_pos = fit_func_inv(pos, *dpd_peak_pos[rho])
#    pred_pos = fit_func(pred_as, *dpd_peak_pos[rho])
except KeyError:
    sys.exit("Coefficients for density %i not found." % rho)

try:
    orig_as = np.array([subst.a[sol] / rho * nm**(2/3) for nm in Nm])
except KeyError:
    sys.exit("Substance '%s' not defined." % sol)

df = pd.DataFrame({"Nm": Nm, "orig_as": orig_as, "pred_as": pred_as})
df = df.sort_values(["Nm"])
print(df)
fname = "new_as.out"
df.to_csv(fname, sep=",")
print("Data frame saved in %s." % fname)

#V0 = 18 * 1.67e-27 / 1000 * 1e30
#rc_DPD = np.array([(rho * nm * V0)**(1/3) for nm in Nm])
#print("MD peaks position:\n", pos / rc_DPD)
#print("DPD peaks position:\n", pred_pos)


