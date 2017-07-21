#!/usr/bin/env python
"""Usage:
    process_surf_ten.py <outfile> [--plot_B_cut --plot_A_cut --plot3d]

Load surface tension, create a data frame and plot results.

Options:
    --plot    Plot cuts holding A or B constant
    --plot3d  Plot gamma vs A and B in 3d

07/06/17
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from docopt import docopt


def plot_B_cut(cut_df, Bval):
    plt.cla()
    plt.plot(cut_df.A, cut_df.gamma, "+-", ms=10)
    plt.xlabel("$A$")
    plt.ylabel("$\gamma$")
    plt.title("$B = %s$" % Bval)
    figname = "cut_B%i.png" % Bval
    plt.savefig(figname, bbox_inches="tight")
    print("Plot saved in %s." % figname)


def plot_A_cut(cut_df, Aval):
    plt.cla()
    plt.plot(cut_df.B, cut_df.gamma, "+-", ms=10)
    plt.xlabel("$B$")
    plt.ylabel("$\gamma$")
    plt.title("$A = %s$" % Aval)
    figname = "cut_A%i.png" % Aval
    plt.savefig(figname, bbox_inches="tight")
    print("Plot saved in %s." % figname)


def plot_3d(df):
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    pass # unfinished


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

#print(df)
dfname = dirn + "/gamma_list_cleared.out"
df = df[["rd", "A", "B", "gamma", "std"]]
df.to_csv(dfname, sep=",")
print("Cleared data frame saved in %s." % dfname)
df2 = df[["A", "B", "gamma"]]
df2.to_csv(dirn + "/for_gnuplot.out", sep="\t", header=False)

As = np.array(list(set(df.A)))
Bs = np.array(list(set(df.B)))

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
       

