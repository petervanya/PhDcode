#!/usr/bin/env python
"""Usage:
    process_densities.py <outfile> [--plot_B_cut --plot_A_cut --plot-imshow]

Load the raw density file, create a data frame out of it and plot results.

Options:
    --plot_B_cut     Plot slices with constant B
    --plot_A_cut     Plot slices with constant A
    --plot-imshow    Plot 3D image

21/07/17
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from docopt import docopt


def plot_B_cut(cut_df, Bval):
    plt.clf()
    plt.plot(cut_df.A, cut_df.rho, "+-", ms=10)
    plt.xlim([-105, 0])
    plt.ylim([0, 13])
    plt.grid()
    plt.xlabel("$A$")
    plt.ylabel("$\\rho$")
    plt.title("$Fixed B$ = %s" % Bval)
    figname = "rho_fixed_B%i.png" % Bval
    plt.savefig(figname, bbox_inches="tight")
    print("Plot saved in %s." % figname)


def plot_A_cut(cut_df, Aval):
    plt.clf()
    c = np.polyfit(np.log(cut_df.B), np.log(cut_df.rho), 1)
    plt.plot(cut_df.B, cut_df.rho, "+-", ms=10)
    plt.xlim([0, 105])
    plt.ylim([0, 13])
    plt.grid()
    plt.xlabel("$B$")
    plt.ylabel("$\\rho$")
    plt.title("Fixed $A$ = %s | log fit: $%.3f x^{%.3f}$" % (Aval, c[1], c[0]))
    figname = "rho_fixed_A%i.png" % Aval
    plt.savefig(figname, bbox_inches="tight")
    print("Plot saved in %s." % figname)

#    plt.clf()
#    print(c)
#    B_log = np.polyval(c, np.log(cut_df.B))
#    plt.plot(np.log(cut_df.B), np.log(cut_df.rho), "+-")
#    plt.plot(np.log(cut_df.B), B_log)
#    plt.title("Fixed $A$ = %s | log fit: $%.3f x^{%.3f}$" % (Aval, c[1], c[0]))
#    figname = "log_rho_cut_A%i.png" % Aval
#    plt.savefig(figname, bbox_inches="tight")
#    print("Log plot saved in %s." % figname)


args = docopt(__doc__)
outfile = args["<outfile>"]
print("Reading %s..." % outfile)
df = pd.read_csv(outfile, sep=" ", header=None)
df.columns = ["sys", "rho"]
df["rd"] = [float(line.split("_")[1]) for line in df.sys]
df["A"] = [float(line.split("_")[3]) for line in df.sys]
df["B"] = [float(line.split("_")[5]) for line in df.sys]
df = df.drop("sys", 1)
df = df[df.A < 0]

df = df[["rd", "A", "B", "rho"]]
dfname = "density_cleared.out"
df.to_csv(dfname)
print("Cleared data frame saved in %s." % dfname)
df2 = df[["A", "B", "rho"]]
df2.to_csv("density_gnuplot.out", sep="\t", header=False)

As = np.array(sorted(list(set(df.A))))
Bs = np.array(sorted(list(set(df.B))))

if args["--plot_B_cut"]:
    for B in Bs:
        cut_df = df[["A", "rho"]][df.B == B]
        cut_df = cut_df.sort_values(["A"])
        plot_B_cut(cut_df, B)

if args["--plot_A_cut"]:
    for A in As:
        cut_df = df[["B", "rho"]][df.A == A]
        cut_df = cut_df.sort_values(["B"])
        plot_A_cut(cut_df, A)

if args["--plot-imshow"]:
    M, N = len(As), len(Bs)
    As_d = {k: v for k, v in zip(sorted(As), range(M))}
    Bs_d = {k: v for k, v in zip(sorted(Bs), range(N))}
    extent = np.array([5, 100, -5, -100]) + np.array([-1, 1, 1, -1]) * 2.5
    rho = np.zeros((M, N))
    for i, r in df2.iterrows():
        rho[As_d[r["A"]], Bs_d[r["B"]]] = r["rho"]
    plt.imshow(rho, extent=extent)
    plt.xlabel("$B$")
    plt.ylabel("$A$")
    plt.colorbar()
    plt.savefig("density_imshow.png")


