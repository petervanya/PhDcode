#!/usr/bin/env python
"""Usage:
    process_densities.py <outfile> [--plot_B_cut --plot_A_cut]

Load the raw density file, create a data frame out of it and plot results.

21/07/17
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
from docopt import docopt


def plot_B_cut(cut_df, Bval):
    plt.cla()
    plt.plot(cut_df.A, cut_df.rho, "+-", ms=10)
    plt.xlim([-105, 0])
    plt.ylim([0, 20])
    plt.grid()
    plt.xlabel("$A$")
    plt.ylabel("$\\rho$")
    plt.title("$B = %s$" % Bval)
    figname = "rho_cut_B%i.png" % Bval
    plt.savefig(figname, bbox_inches="tight")
    print("Plot saved in %s." % figname)


def plot_A_cut(cut_df, Aval):
    plt.cla()
    plt.plot(cut_df.B, cut_df.rho, "+-", ms=10)
    plt.xlim([0, 105])
    plt.ylim([0, 20])
    plt.grid()
    plt.xlabel("$B$")
    plt.ylabel("$\\rho$")
    plt.title("$A = %s$" % Aval)
    figname = "rho_cut_A%i.png" % Aval
    plt.savefig(figname, bbox_inches="tight")
    print("Plot saved in %s." % figname)


args = docopt(__doc__)
outfile = args["<outfile>"]
print("Reading %s..." % outfile)
df = pd.read_csv(outfile, sep=" ", header=None)
df.columns = ["sys", "rho"]
df["rd"] = [float(line.split("_")[1]) for line in df.sys]
df["A"] = [float(line.split("_")[3]) for line in df.sys]
df["B"] = [float(line.split("_")[5]) for line in df.sys]
df = df.drop("sys", 1)

df = df[["rd", "A", "B", "rho"]]
dfname = "density_cleared.out"
df.to_csv(dfname, sep=",")
print("Cleared data frame saved in %s." % dfname)
df2 = df[["A", "B", "rho"]]
df2.to_csv("density_gnuplot.out", sep="\t", header=False)

As = np.array(list(set(df.A)))
Bs = np.array(list(set(df.B)))

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


