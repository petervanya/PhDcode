#!/usr/bin/env python
"""Usage:
    fit_gamma_A_B.py <frame> <func> [--iter <it> --plot_A_cut --plot_B_cut]

Fit surface tension as a function of A and B via curve_fit.

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
ymax = {0.55: 150, 0.65: 120, 0.75: 120, 0.85: 120}


def func1(X, c1, c2, c3, c4, c5):
    """(c1 * A**2 + c2 * A + c3) * (B - c4)**c5"""
    A, B = X
    return (c1 * A**2 + c2 * A + c3) * (B - c4)**c5


def func2(X, c1, c2, c3, c4, c5, c6):
    """(c1 * A**2 + c2 * A + c3) * (B - c4 + c5 * A)**c6"""
    A, B = X
    return (c1 * A**2 + c2 * A + c3) * (B - c4 + c5 * A)**c6


def func3(X, c1, c2, c3, c4, c5, c6):
    """(c1 * A**2 + c2 * A + c3) * (B - c4)**(c5 + c6 * A)"""
    A, B = X
    return (c1 * A**2 + c2 * A + c3) * (B - c4)**(c5 + c6 * A)


def func4(X, c1, c2, c3, c4):
    """(c1 * A**2 + c2 * A) * (B - c3)**c4"""
    A, B = X
    return (c1 * A**2 + c2 * A) * (B - c3)**c4


def func5(X, c1, c2, c3, c4):
    """(c1 * A**2 + c2) * (B - c3)**c4"""
    A, B = X
    return (c1 * A**2 + c2) * (B - c3)**c4


def func6(X, c1, c2, c3, c4, c5, c6):
    """(c1 * A**3 + c2 * A**2 + c3 * A + c4) * (B - c5)**c6"""
    A, B = X
    return (c1 * A**3 + c2 * A**2 + c3 * A + c4) * (B - c5)**c6


def func7(X, c1, c2, c3, c4, c5):
    """(c1 * A**3 + c2 * A**2 + c3 * A) * (B - c4)**c5"""
    A, B = X
    return (c1 * A**3 + c2 * A**2 + c3 * A) * (B - c4)**c5


def func8(X, c1, c2, c3, c4):
    """(c1 * A**2 + c2 * A) * exp(-c3 * B + c4)"""
    A, B = X
    return (c1 * A**2 + c2 * A) * np.exp(-c3 * B + c4)


def plot_B_cut(cut_df, Bval, popt):
    plt.clf()
    gammafit = np.array([func((xi, Bval), *popt) for xi in cut_df.A])
    plt.plot(cut_df.A, cut_df.gamma, "+-", ms=10)
    plt.plot(cut_df.A, gammafit, ms=10)
    plt.xlim([-105, 0])
    plt.ylim([-1, ymax[rd]])
    plt.grid()
    plt.xlabel("$A$")
    plt.ylabel("$\\gamma$")
    plt.title("Fixed $B$ = %s" % Bval)
    figname = "fit_gamma_fixed_B%i.png" % Bval
    plt.savefig(figname, bbox_inches="tight")
    print("Plot saved in %s." % figname)


def plot_A_cut(cut_df, Aval, popt):
    plt.clf()
    gammafit = np.array([func((Aval, xi), *popt) for xi in cut_df.B])
    plt.plot(cut_df.B, cut_df.gamma, "+-", ms=10)
    plt.plot(cut_df.B, gammafit, ms=10)
    plt.xlim([0, 105])
    plt.ylim([-1, ymax[rd]])
    plt.grid()
    plt.xlabel("$B$")
    plt.ylabel("$\\gamma$")
    plt.title("Fixed $A$ = %s" % Aval)
    figname = "fit_gamma_fixed_A%i.png" % Aval
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

df.columns = ["sys", "gamma", "dg"]
df["rd"] = [float(line.split("_")[1]) for line in df.sys]
df["A"] = [float(line.split("_")[3]) for line in df.sys]
df["B"] = [float(line.split("_")[5]) for line in df.sys]
rd = df.rd.values[0]
df = df[df.B > Bcut[rd]]  # do not take too small B values
XA = df.A.values
XB = df.B.values
gamma = df.gamma.values
f = int(args["<func>"])
Nf = 8
if f not in range(1, Nf+1):
    sys.exit("Choose function <f> from %s." % range(1, Nf+1))

print("===== Fitting MDPD surface tension w.r.t. A, B =====")

if f == 1:
    func = func1
    guess = [0.06, -0.2, -10.0, 5.0, -0.8]
elif f == 2:
    func = func2
    guess = [0.06, -0.2, -10.0, 5.0, 0.001, -0.8]
elif f == 3:
    func = func3
    guess = [0.06, -0.2, -10.0, 5.0, -0.8, 0.1]
elif f == 4:
    func = func4
    guess = [0.06, -0.2, 5.0, -0.8]
elif f == 5:
    func = func5
    guess = [0.06, -0.2, -10.0, -0.8]
if f == 6:
    func = func6
    guess = [0.01, 0.06, -0.2, -10.0, 5.0, -0.8]
if f == 7:
    func = func7
    guess = [0.01, -0.2, -10.0, 5.0, -0.8]
if f == 8:
    func = func8
    guess = [0.06, -0.2, 1.0, 1.0]

print(func.__doc__)
popt, pcov = curve_fit(func, (XA, XB), gamma, p0=guess, maxfev=Nit)
print("Guess, params, errors, ratios")
np.set_printoptions(precision=6, suppress=True)
print(np.c_[guess, popt, np.diag(pcov), np.abs(np.diag(pcov) / popt)])

As = np.array(sorted(list(set(df.A))))
Bs = np.array(sorted(list(set(df.B))))

if args["--plot_B_cut"]:
    for B in Bs:
        cut_df = df[["A", "gamma"]][df.B == B]
        cut_df = cut_df.sort_values(["A"])
        plot_B_cut(cut_df, B, popt)

if args["--plot_A_cut"]:
    for A in As:
        cut_df = df[["B", "gamma"]][df.A == A]
        cut_df = cut_df.sort_values(["B"])
        plot_A_cut(cut_df, A, popt)


