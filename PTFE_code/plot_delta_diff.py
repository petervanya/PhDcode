#!/usr/bin/env python
"""Usage:
    plot_delta_diff.py [--tex --fmt <fmt>]

Plot difference between parallel and normal diffusivity

Options:
    --tex          Use tex for axes
    --fmt <fmt>    Plot format [default: png]

09/05/17, modified 17/10/17
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import glob, sys
from docopt import docopt

matplotlib.rcParams.update({'font.size': 24})


def plot_el(df, el, widths, lmbdas, figfmt):
    """Plot parallel and normal diffusivities for electrode
    * df: dataframe
    * el: electrode
    * widths: range of widths
    * lmbdas: range of lmbdas
    """
    fig, ax = plt.subplots(figsize=(8, 6))
    linestyles = ["-","--",":","-."]

    for (d, ls) in zip(widths, linestyles):
        lbl = "%i nm" % d
        lbl = "%i nm" % d
        Dx = np.array(df["Dx"][(df.d == d) & (df.elect == el)])
        Dx_err = np.array(df["Dx_err"][(df.d == d) & (df.elect == el)])
        Dyz = np.array(df["Dyz"][(df.d == d) & (df.elect == el)])
        Dyz_err = np.array(df["Dyz_err"][(df.d == d) & (df.elect == el)])
        D3d = np.array(df["D3d"][(df.d == d) & (df.elect == el)])
        D3d_err = np.array(df["D3d_err"][(df.d == d) & (df.elect == el)])
        eb = Dx_err / Dyz + Dx / Dyz**2 * Dyz_err   # error for (Dx - Dyz) / Dyz
#        eb = Dyz_err / Dx + Dyz / Dx**2 * Dx_err    # error for (Dx - Dyz) / Dx
#        eb = Dyz_err / D3d + Dyz / D3d**2 * D3d_err \
#                + Dx_err / D3d + Dx / D3d**2 * D3d_err   # error for (Dx - Dyz) / D3d
#        plt.errorbar(lmbdas, (Dyz - Dx) / Dx, yerr=eb, \
        plt.errorbar(lmbdas, (Dyz - Dx) / Dyz, yerr=eb, \
                fmt=ls, lw=2, ms=2, mew=2, capsize=3, label=lbl)

    plt.xlim([lmbdas[0]-1, lmbdas[-1]+1])
    plt.xticks([4, 8, 12, 16, 20, 24])
    plt.ylim([0, 0.5])
    plt.xlabel("$\lambda$")
#    plt.ylabel("($D_{\mathrm{yz}} - D_{\mathrm{x}}) / D_{\mathrm{x}}$")
    plt.ylabel("($D_{\|} - D_{\\bot}) / D_{\|}$")

    plt.grid()
    plt.legend(loc="best")
    plt.title("Anisotropy, " + el.lower())
    plotname = "delta_%s.%s" % (el.lower(), figfmt)
    plt.savefig(plotname, bbox_inches='tight')


if __name__ == "__main__":
    args = docopt(__doc__)
    elects = ["Carbon", "Quartz"] #
    lmbdas = list(range(4, 25, 2))
    widths = [5, 10, 15, 20]
    print("Water uptakes:", lmbdas)
    print("Widths:", widths)
    figfmt = args["--fmt"]
    if args["--tex"]:
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
    
    Nd = len(lmbdas) * len(elects) * len(widths)
    df_name = "diff_avg.csv"
    try:
        master_df = pd.read_csv(df_name).dropna()
    except FileNotFoundError:
        sys.exit("Dataframe %s not found." % df_name)
   
    for el in elects:
        plot_el(master_df, el, widths, lmbdas, figfmt)


