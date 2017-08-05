#!/usr/bin/env python
"""Usage:
    fit_density.py <profile> [--cut <rhoc> --plot --saveplot]

Fit the density profile using hyperbolic tan.

Options:
    --cut <rhoc>   Cutoff below which density is zero [default: 0.5]

20/07/17
"""
import numpy as np
from scipy.optimize import curve_fit
import sys
from docopt import docopt


def func(x, c1, c2, c3, c4, c5):
    return c1 * (-np.tanh(c2 * (np.abs(x - c3) - c4)) + 1.0) / 2.0 + c5


args = docopt(__doc__)
cut = float(args["--cut"])
pr = args["<profile>"]
try:
    A = np.loadtxt(pr)
except FileNotFoundError:
    sys.exit("File %s not found." % pr)

x, rho = A[:, 0], A[:, 1]
dx = x[1] - x[0]
N = len(x)
L = N * dx
guess = [5.0, 10.0, L / 2, 2.0, 0.0]
popt, pcov = curve_fit(func, x, rho, p0=guess)

print("===== Fitting the MDPD density profile with tanh =====")
print("Params: %s" % popt)
print("Fitted density: %.3f" % (popt[0] + popt[4]))
print("Interface slope: %.3f" % popt[1])
print("Interface width: %.3f" % (1.0 / popt[1]))

if args["--plot"]:
    import matplotlib.pyplot as plt
    rhofit = np.array([func(xi, *popt) for xi in x])
    plt.plot(x, rho, label="profile")
    plt.plot(x, rhofit, label="fit")
    plt.xlabel("$x$")
    plt.ylabel("$\\rho$")
    plt.xlim([0, L])
    plt.ylim(bottom=0.0)
    plt.legend()
    if args["--saveplot"]:
        figname = "density_fit.png"
        plt.savefig(figname)
        print("Plot saved in %s." % figname)
    else:
        plt.show()


