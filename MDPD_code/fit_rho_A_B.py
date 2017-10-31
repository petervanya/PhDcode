#!/usr/bin/env python
"""Usage:
    fit_rho_A_B.py <frame> [--iter <it> --plot_A_cut --plot_B_cut]

Fit density as a function of A and B via curve_fit.

Options:
    --iter <it>    Number of iterations [default: 10000]
    --plot         Plot the curves

23/10/17
"""
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import sys
from docopt import docopt


Bcut = {0.55: 20, 0.65: 15, 0.75: 10, 0.85: 10}
ymax = {0.55: 30, 0.65: 20, 0.75: 12, 0.85: 12}


def func(X, c1, c2, c3, c4):
    A, B = X
    return c1 + c2 * (-A) * (B - c3)**c4


def plot_B_cut(cut_df, Bval, popt):
    plt.clf()
    rhofit = np.array([func((xi, Bval), *popt) for xi in cut_df.A])
    plt.plot(cut_df.A, cut_df.rho, "+-", ms=10)
    plt.plot(cut_df.A, rhofit, ms=10)
    plt.xlim([-105, 0])
    plt.ylim([0, ymax[rd]])
    plt.grid()
    plt.xlabel("$A$")
    plt.ylabel("$\\rho$")
    plt.title("Fixed $B$ = %s" % Bval)
    figname = "fit_rho_fixed_B%i.png" % Bval
    plt.savefig(figname, bbox_inches="tight")
    print("Plot saved in %s." % figname)


def plot_A_cut(cut_df, Aval, popt):
    plt.clf()
    rhofit = np.array([func((Aval, xi), *popt) for xi in cut_df.B.values])
    plt.plot(cut_df.B, cut_df.rho, "+-", ms=10)
    plt.plot(cut_df.B, rhofit, ms=10)
    plt.xlim([0, 105])
    plt.ylim([0, ymax[rd]])
    plt.grid()
    plt.xlabel("$B$")
    plt.ylabel("$\\rho$")
    plt.title("Fixed $A$ = %s" % Aval)
    figname = "fit_rho_fixed_A%i.png" % Aval
    plt.savefig(figname, bbox_inches="tight")
    print("Plot saved in %s." % figname)


args = docopt(__doc__)
fname = args["<frame>"]
Nit = int(args["--iter"])
print("Reading %s..." % fname)
try:
    df = pd.read_csv(fname, sep=" ", header=None)
except FileNotFoundError:
    sys.exit("File %s not found." % fname)

df.columns = ["sys", "rho"]
df["rd"] = [float(line.split("_")[1]) for line in df.sys]
df["A"] = [float(line.split("_")[3]) for line in df.sys]
df["B"] = [float(line.split("_")[5]) for line in df.sys]
rd = df.rd.values[0]
df = df[df.B > Bcut[rd]]  # do not take too small B values
XA = df.A.values
XB = df.B.values
rho = df.rho.values

print("===== Fitting MDPD density w.r.t. A, B =====")
print("rho(A,B) = c1 + c2 (-A) * (B - c3)**c4")
guess = [3.0, 1.0, 1.0, -0.8]
popt, pcov = curve_fit(func, (XA, XB), rho, p0=guess, maxfev=Nit)
print("Guess, params, errors, ratios")
np.set_printoptions(precision=6, suppress=True)
print(np.c_[guess, popt, np.diag(pcov), np.abs(np.diag(pcov) / popt)])

As = np.array(sorted(list(set(df.A))))
Bs = np.array(sorted(list(set(df.B))))

if args["--plot_B_cut"]:
    for B in Bs:
        cut_df = df[["A", "rho"]][df.B == B]
        cut_df = cut_df.sort_values(["A"])
        plot_B_cut(cut_df, B, popt)

if args["--plot_A_cut"]:
    for A in As:
        cut_df = df[["B", "rho"]][df.A == A]
        cut_df = cut_df.sort_values(["B"])
        plot_A_cut(cut_df, A, popt)


