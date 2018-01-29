#!/usr/bin/env python
"""Usage: 
    plot_avg_diff.py <methods> [--tex --fmt <fmt>]

Plot diffusivity from all seeds, e.g.
$ ./plot_avg_diff.py "../Method_3_s123*"

Options:
    --tex          Use tex for axes
    --fmt <fmt>    Plot format [default: png]

05/05/17
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import glob, sys
from docopt import docopt

matplotlib.rcParams.update({'font.size': 36})
#D0 = 1.752226e-08
#D0 = 2.068688
D0 = 0.3


def plot_el(df, el, D, widths, lmbdas, figfmt):
    """Plot parallel and normal diffusivities for electrode
    * df: dataframe
    * el: electrode
    * D: either 'Dx' or 'Dyz'
    * widths: range of widths
    * lmbdas: range of lmbdas
    """
    fig, ax = plt.subplots(figsize=(8, 6))
    linestyles = ["-","--",":","-."]

    for (d, ls) in zip(widths, linestyles):
        lbl = "%i nm" % d
        sel = np.array(df[["lmbda", D]][(df.d == d) & (df.elect == el)])
        eb = np.array(df[D+"_err"][(df.d == d) & (df.elect == el)])
        plt.errorbar(sel[:, 0], sel[:, 1] / D0, yerr=eb / D0, \
                fmt=ls, lw=2, ms=2, mew=2, capsize=3.0, label=lbl)

    plt.xlim([lmbdas[0]-1, lmbdas[-1]+1])
    plt.xticks([4, 8, 12, 16, 20, 24])
    plt.ylim([0, 1.0]) #0.2])
    plt.yticks(np.linspace(0, 1, 6))
#    plt.yticks([0, 0.05, 0.10, 0.15, 0.20])
    plt.xlabel("$\lambda$")
#    ylbl = {"Dx": "Normal", "Dyz": "Parallel", "D3d": "3d"}
#    plt.ylabel("$D_{\mathrm{%s}} \; (10^{-9} \; \mathrm{m^2/s}) $" % ylbl[D])
    ylbl = {"Dx": "\\bot", "Dyz": "\|", "D3d": "3D"}
    plt.ylabel("$D_{\mathrm{%s}} \,/\, D_0 $" % ylbl[D])

    plt.grid()
    plt.legend(loc="best", prop={'size': 24})
    ttl = {"Dx": "Normal", "Dyz": "Parallel", "D3d": "3D"}
    plt.title(ttl[D] + ", " + el.lower())
    plotname = "%s_avg_%s.%s" % (D, el.lower(), figfmt)
    plt.savefig(plotname, bbox_inches='tight')


if __name__ == "__main__":
    args = docopt(__doc__)
    dirs = glob.glob(args["<methods>"])
    print(dirs)
    elects = ["Carbon", "Quartz"] #
    lmbdas = list(range(4, 25, 2))
    widths = [5, 10, 15, 20]
    print("Water uptakes:", lmbdas)
    print("Widths:", widths)
    figfmt = args["--fmt"]
    if args["--tex"]:
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
    
    Nd = len(lmbdas) * len(widths) * len(elects)
    Dx_mat = np.zeros((Nd, len(dirs)))
    Dyz_mat = np.zeros((Nd, len(dirs)))
    D3d_mat = np.zeros((Nd, len(dirs)))
    
    for i in range(len(dirs)):
        df = pd.read_csv(dirs[i] + "/DataCommon/all_diff.csv").dropna()
        if len(df) != Nd:
            print("Length of data frame should be equal to %i." % Nd)
        Dx_mat[:, i] = df.Dx
        Dyz_mat[:, i] = df.Dyz
        D3d_mat[:, i] = df.D3d
    
    Dx_avg = np.average(Dx_mat, 1)
    Dx_err = np.std(Dx_mat, 1)
    Dyz_avg = np.average(Dyz_mat, 1)
    Dyz_err = np.std(Dyz_mat, 1)
    D3d_avg = np.average(D3d_mat, 1)
    D3d_err = np.std(D3d_mat, 1)
    
    master_df = pd.read_csv(dirs[0] + "/DataCommon/all_diff.csv").dropna()
    master_df = master_df[["elect", "d", "lmbda"]]
    master_df["Dx"] = Dx_avg
    master_df["Dx_err"] = Dx_err
    master_df["Dyz"] = Dyz_avg
    master_df["Dyz_err"] = Dyz_err
    master_df["D3d"] = D3d_avg
    master_df["D3d_err"] = D3d_err
    master_df.to_csv("diff_avg.csv")
    print("Data frame saved in diff_avg.csv.")
   
    for el in elects:
        plot_el(master_df, el, "Dx", widths, lmbdas, figfmt)
        plot_el(master_df, el, "Dyz", widths, lmbdas, figfmt)
        plot_el(master_df, el, "D3d", widths, lmbdas, figfmt)


