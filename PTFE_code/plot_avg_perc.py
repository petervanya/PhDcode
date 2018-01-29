#!/usr/bin/env python
"""Usage:
    plot_avg_perc.py <methods> [--tex --fmt <fmt>]

Plot diffusivity from all seeds, e.g.
$ ./plot_avg_perc.py "../Method_3_s123*"

Options:
    --tex          Use tex for axes
    --fmt <fmt>    Plot format [default: png]

05/05/17
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import glob, sys
from docopt import docopt
import pandas as pd

matplotlib.rcParams.update({'font.size': 36})
r_DPD = 8.14e-10


def plot_el(df, el, dim, widths, lmbdas, figfmt):
    """Plot parallel and normal diffusivities for electrode
    * df: dataframe
    * el: electrode
    * dim: '2d' or '3d'
    * widths: range of widths
    * lmbdas: range of lmbdas
    """
    fig, ax = plt.subplots(figsize=(8, 6))
    Pinf = "P" + dim
    linestyles = ["-","--",":","-."]

    for (d, ls) in zip(widths, linestyles):
        lbl = "%i nm" % d
        sel = np.array(df[["lmbda", Pinf]][(df.d == d) & (df.elect == el)])
        eb = np.array(df[Pinf + "_err"][(df.d == d) & (df.elect == el)])
        if Pinf == "P3d":
            scale = (d + 10.) / (40. * 0.814)  # rescale w.r.t. amount of water
            plt.errorbar(sel[:, 0], sel[:, 1] / scale, yerr=eb, \
                    fmt=ls, lw=2, ms=2, mew=2, capsize=3, label=lbl)
        else:
            plt.errorbar(sel[:, 0], sel[:, 1], yerr=eb, \
                fmt=ls, lw=2, ms=2, mew=2, capsize=3, label=lbl)

    plt.xlim([lmbdas[0]-1, lmbdas[-1]+1])
    plt.xticks([4, 8, 12, 16, 20, 24])
    plt.ylim([0, 1])
    plt.xlabel("$\lambda$")
    plt.ylabel("$P_{\infty}$")

    plt.grid()
    if dim == "2d":
        plt.legend(loc=4, prop={'size': 24})
    else:
        plt.legend(loc=2, prop={'size': 24})
#    plt.title(dim.upper() + " percolation, thin film, " + el.lower())
    plt.title(dim.upper() + " percolation, " + el.lower())
    plotname = "Pinf_%s_avg_%s.%s" % (dim, el.lower(), figfmt)
    plt.savefig(plotname, bbox_inches='tight')


if __name__ == "__main__":
    args = docopt(__doc__)
    dirs = glob.glob(args["<methods>"])
    print(dirs)
    elects = ["Carbon", "Quartz"]
    lmbdas = list(range(4, 25, 2))
    widths = [5, 10, 15, 20]
    print("Water uptakes:", lmbdas)
    print("Widths:", widths)
    figfmt = args["--fmt"]
    if args["--tex"]:
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
    
    Nd = len(lmbdas) * len(elects) * len(widths)
    P2d_mat = np.zeros((Nd, len(dirs)))
    P3d_mat = np.zeros((Nd, len(dirs)))
    
    for i in range(len(dirs)):
        df = pd.read_csv(dirs[i] + "/DataCommon/all_perc.csv").dropna()
        P2d_mat[:, i] = df.P2d
        P3d_mat[:, i] = df.P3d
    
    P2d_avg = np.average(P2d_mat, 1)
    P2d_err = np.std(P2d_mat, 1)
    P3d_avg = np.average(P3d_mat, 1)
    P3d_err = np.std(P3d_mat, 1)

    master_df = pd.read_csv(dirs[0] + "/DataCommon/all_perc.csv").dropna()
    master_df = master_df[["elect", "d", "lmbda"]]
    master_df["P2d"] = P2d_avg
    master_df["P2d_err"] = P2d_err
    master_df["P3d"] = P3d_avg
    master_df["P3d_err"] = P3d_err
    master_df.to_csv("perc_avg.csv")
    print("Data frame saved in perc_avg.csv.")

    for el in elects:
        plot_el(master_df, el, "2d", widths, lmbdas, figfmt)
        plot_el(master_df, el, "3d", widths, lmbdas, figfmt)

