#!/usr/bin/env python
"""Usage:
    fit_dpd_rdfs.py <rho> <a> <fitfunc> [--bounds --it <it>]

Fit RDF by one function. Fitting functions:
    * cv (Einstein heat cap rectification)
    * exp (simple exponential)
    * expsq (gaussian)
    * expsq_6 (gaussian with additional parameter)

Options:
    --it <it>      Number of iterations [default: 20000]
    --bounds       Include bounds on fitting parameters

pv278@cam.ac.uk, 08/06/17
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sys, os
from docopt import docopt


def func_cv(x, a1, a2, a3, a4, a5):
    """'Basis' functions:
    * cv (Einstein heat cap term)
    * exp(-x**2)
    * sin(x)
    Result: not good"""
    cv = (a5 / x)**2 * np.exp(a5 / x) / (np.exp(a5 / x) - 1.0)**2
    f = a1 * x**3 * np.exp(-a2 * x**2) * (np.sin(a3 * x + a4) + 1.0) * cv
    return f


def func_exp(x, a1, a2, a3, a4, a5):
    """'Basis' functions: exp(-abs(x)), sin(x), tanh(x)
    5 parameters.
    Result: good, but main peak not captured and wrong frequency"""
    rect = (np.tanh(a4 * (x - a5)) + 1.0) / 2
    f = (a1 * np.exp(-np.abs(x - a2)) * np.sin(a2 * (x - a3)) + 1.0) * rect
    return f


def func_expsq(x, a1, a2, a3, a4, a5):
    """Functions: exp(-x**2), sin(x), tanh(x)
    5 parameters
    Result: sometimes good, sometimes very bad"""
    rect = (np.tanh(a4 * (x - a5)) + 1.0) / 2
    f = (a1 * np.exp(-(x - a2)**2) * np.sin(a3 * (x - a2)) + 1.0) * rect
    return f


def func_expsq_relax(x, a1, a2, a3, a4, a5, a6):
    """Functions: exp(-x**2), sin(x), tanh(x)
    6 parameters, relax exp and sin shifts
    Result: sometimes good, sometimes very bad"""
    rect = (np.tanh(a4 * (x - a5)) + 1.0) / 2
    f = (a1 * np.exp(-(x - a2)**2) * np.sin(a3 * (x - a6)) + 1.0) * rect
    return f


def func_expsq_6(x, a1, a2, a3, a4, a5, a6):
    """Functions: exp(-x**2), sin(x), tanh(x)
    6 parameters
    Result: not much better from expsq with 5 parameters"""
    rect = (np.tanh(a4 * (x - a5)) + 1.0) / 2
    f = (a1 * np.exp(- a6 * (x - a2)**2) * np.sin(a3 * (x - a2)) + 1.0) * rect
    return f


def func_expsq_6_relax(x, a1, a2, a3, a4, a5, a6, a7):
    """Functions: exp(-x**2), sin(x), tanh(x)
    7 parameters, relax exp and sin shifts
    Result: not much better from expsq with 5 parameters"""
    rect = (np.tanh(a4 * (x - a5)) + 1.0) / 2
    f = (a1 * np.exp(- a6 * (x - a2)**2) * np.sin(a3 * (x - a7)) + 1.0) * rect
    return f


def plot_fit(r, rdf, yfit, rho, a, fitfunc, bounds):
    plt.clf()
    plt.plot(r, rdf, lw=3)
    plt.plot(r, yfit, lw=1.5)
    plt.title("$\\rho = %i, a = %i$, %s" % (rho, a, fitfunc))
    plt.xlabel("$r$ (DPD units)")
    plt.ylabel("$g(r)$")
    plt.xlim([0, 5])
    plt.ylim([0, 2])
    plt.grid()
    plotname = "rdf_fit_rho_%i_a_%i.png" % (rho, a)
    if bounds:
        plotname = "rdf_fit_rho_%i_a_%i_b.png" % (rho, a)
    savedir = "Fit_" + fitfunc
    path = savedir + "/" + plotname
    plt.savefig(path)
    print(path + " fit plotted.")


if __name__ == "__main__":
    args = docopt(__doc__)
    rho = float(args["<rho>"])
    a = float(args["<a>"])
    fitfunc = args["<fitfunc>"]
    
    fname = "rho_%i_a_%i" % (rho, a) + "/rdf_1.out"
    func_dir = {"cv": func_cv, "exp": func_exp, "expsq": func_expsq, \
            "expsq_6": func_expsq_6, "expsq_relax": func_expsq_relax,
            "expsq_6_relax": func_expsq_6_relax}
    if fitfunc not in func_dir.keys():
        sys.exit("Choose fitting function from the list:\n%s" \
                % list(func_dir.keys()))

    savedir = "Fit_" + fitfunc
    if not os.path.isdir(savedir):
        os.makedirs(savedir)

    print("===== Fitting RDF %s" % "with bounds" if args["--bounds"] else "")
    print("rho: %i | a: %i | function: %s" % (rho, a, fitfunc))
    
    try:
        A = np.loadtxt(fname)
    except FileNotFoundError:
        print("%s not found." % fname)
    
    r, rdf = A[:, 0], A[:, 1]
    guess = [0.5, 0.85, 10.0, 1.0, 2.0]
    mins = [-5.0, -5.0, 0.0,  0.0, 0.0]
    maxs = [5.0, 10.0, 20.0, 20.0, 1.0]
    if fitfunc in ["expsq_6", "expsq_relax"]:
        guess = [0.5, 0.85, 10.0, 1.0, 2.0, 1.0]
        mins = [-5.0, -5.0, 0.0,  0.0, 0.0, 0.0]
        maxs = [5.0, 10.0, 20.0, 20.0, 2.0, 3.0]
    if fitfunc in ["expsq_6_relax"]:
        guess = [0.5, 0.85, 10.0, 1.0, 2.0, 1.0, 1.0]

    
    Nit = int(args["--it"])
    popt, pcov = curve_fit(func_dir[fitfunc], r, rdf, p0=guess, maxfev=Nit)
    if args["--bounds"]:
        popt, pcov = curve_fit(func_dir[fitfunc], r, rdf, p0=guess, \
                maxfev=Nit, bounds=(mins, maxs))
    print("Guess  :", "".join(["%10.6f" % i for i in guess]))
    print("Params :", "".join(["%10.6f" % i for i in popt]))
    yfit = func_dir[fitfunc](r, *popt)
    print("Correlation:", np.corrcoef(rdf, yfit)[0, 1])
    plot_fit(r, rdf, yfit, rho, a, fitfunc, args["--bounds"])
        

