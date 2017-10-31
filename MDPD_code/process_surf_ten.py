#!/usr/bin/env python
"""Usage:
    process_surf_ten.py <outfile> [--plot_B_cut --plot_A_cut --plot-imshow]

Load the raw surface tension, create a data frame and plot results.

Options:
    --plot_B_cut     Plot slices with constant B
    --plot_A_cut     Plot slices with constant A
    --plot-imshow    Plot 3D image

07/06/17
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from docopt import docopt


def plot_B_cut(cut_df, Bval):
    plt.cla()
    plt.plot(cut_df.A, cut_df.gamma, "+-", ms=10)
    plt.xlim([-105, 0])
    plt.ylim([0, 120])
    plt.xlabel("$A$")
    plt.ylabel("$\gamma$")
    plt.grid()
    plt.title("$B = %s$" % Bval)
    figname = "cut_B%i.png" % Bval
    plt.savefig(figname, bbox_inches="tight")
    print("Plot saved in %s." % figname)

#    plt.cla()
#    dff = cut_df[cut_df.gamma > 0.0]
#    popt = np.polyfit(np.log(-dff.A), np.log(dff.gamma), 1)
#    A_log = np.polyval(popt, np.log(-cut_df.A))
#    plt.plot(np.log(-cut_df.A), np.log(cut_df.gamma), "+-")
#    plt.plot(np.log(-cut_df.A), A_log)
#    plt.title("Fixed $B$ = %s | log fit: $%.3f x^{%.3f}$" % \
#            (Bval, popt[1], popt[0]))
#    plt.ylim([-5.0, 5.0])
#    plt.grid()
#    figname = "log_gamma_cut_B%i.png" % Bval
#    plt.savefig(figname, bbox_inches="tight")
#    print("Log plot saved in %s." % figname)


def plot_A_cut(cut_df, Aval):
    plt.cla()
    plt.plot(cut_df.B, cut_df.gamma, "+-", ms=10)
    plt.xlim([0, 105])
    plt.ylim([0, 120])
    plt.xlabel("$B$")
    plt.ylabel("$\gamma$")
    plt.title("$A = %s$" % Aval)
    plt.grid()
    figname = "cut_A%i.png" % Aval
    plt.savefig(figname, bbox_inches="tight")
    print("Plot saved in %s." % figname)

#    plt.cla()
#    dff = cut_df[cut_df.gamma > 0.0]
#    popt = np.polyfit(np.log(dff.B), np.log(dff.gamma), 1)
#    B_log = np.polyval(popt, np.log(cut_df.B))
#    plt.plot(np.log(cut_df.B), np.log(cut_df.gamma), "+-")
#    plt.plot(np.log(cut_df.B), B_log)
#    plt.title("Fixed $A$ = %s | log fit: $%.3f x^{%.3f}$" % \
#            (Aval, popt[1], popt[0]))
#    plt.ylim([-5.0, 5.0])
#    plt.grid()
#    figname = "log_gamma_cut_A%i.png" % Aval
#    plt.savefig(figname, bbox_inches="tight")
#    print("Log plot saved in %s." % figname)


args = docopt(__doc__)
outfile = args["<outfile>"]
dirn = "/".join(outfile.split("/")[0:-1])
print("Reading %s..." % outfile)
df = pd.read_csv(outfile, sep=" ", header=None)
df.columns = ["sys", "gamma", "std"]
df["rd"] = [float(line.split("_")[1]) for line in df.sys]
df["A"] = [float(line.split("_")[3]) for line in df.sys]
df["B"] = [float(line.split("_")[5]) for line in df.sys]
df = df.drop("sys", 1)
df = df[df.A < 0]

#print(df)
df = df[["rd", "A", "B", "gamma", "std"]]
dfname = "gamma_cleared.out"
df.to_csv(dfname)
print("Cleared data frame saved in %s." % dfname)
df2 = df[["A", "B", "gamma"]]
df2.to_csv("gamma_gnuplot.out", sep="\t", header=False)

As = np.array(sorted(list(set(df.A))))
Bs = np.array(sorted(list(set(df.B))))

if args["--plot_B_cut"]:
    for B in Bs:
        cut_df = df[["A", "gamma"]][df.B == B]
        cut_df = cut_df.sort_values(["A"])
        plot_B_cut(cut_df, B)

if args["--plot_A_cut"]:
    for A in As:
        cut_df = df[["B", "gamma"]][df.A == A]
        cut_df = cut_df.sort_values(["B"])
        plot_A_cut(cut_df, A)

if args["--plot-imshow"]:
    M, N = len(As), len(Bs)
    As_d = {k: v for k, v in zip(sorted(As), range(M))}
    Bs_d = {k: v for k, v in zip(sorted(Bs), range(N))}
    extent = np.array([5, 100, -5, -100]) + np.array([-1, 1, 1, -1]) * 2.5
    gamma = np.zeros((M, N))
    for i, r in df2.iterrows():
        gamma[As_d[r["A"]], Bs_d[r["B"]]] = r["gamma"]
    plt.imshow(gamma, extent=extent)
    plt.xlabel("$B$")
    plt.ylabel("$A$")
    plt.colorbar()
    figname = "gamma_imshow.png"
    plt.savefig(figname)
    print("Plt imshow saved in %s." % figname)


