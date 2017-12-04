#!/usr/bin/env python
"""Usage:
    process_mdpd_diff.py <outfile> [--plot_B_cut --plot_A_cut --dcut <dc>]
                                   [--plot-imshow --interp <in>]

Load the raw density file, create a data frame out of it and plot results.

Options:
    --plot_B_cut     Plot slices with constant B
    --plot_A_cut     Plot slices with constant A
    --plot-imshow    Plot 3D image
    --dcut <dc>      Cutoff for diffusivity of a solid [default: 0.0025]
    --interp <in>    Imshow interpolation [default: none]

08/09/17, modified 23/10/17
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
from docopt import docopt

# maximum reasonable diffusivities for a given mb cutoff
top_D3d = {0.55: 0.35, 0.65: 0.4, 0.75: 0.5, 0.85: 0.6}


def plot_B_cut(cut_df, Bval, rd):
    plt.clf()
    plt.plot(cut_df.A, cut_df.D3d, "+-", ms=10)
    plt.xlim([-105, 0])
    plt.ylim([0, top_D3d[rd]])
    plt.grid()
    plt.xlabel("$A$")
    plt.ylabel("$D_3d$")
    plt.title("Fixed $B$ = %s" % Bval)
    figname = "D3d_fixed_B%i.png" % Bval
    plt.savefig(figname, bbox_inches="tight")
    print("Plot saved in %s." % figname)


def plot_A_cut(cut_df, Aval, rd):
    plt.clf()
    plt.plot(cut_df.B, cut_df.D3d, "+-", ms=10)
    plt.xlim([0, 105])
    plt.ylim([0, top_D3d[rd]])
    plt.grid()
    plt.xlabel("$B$")
    plt.ylabel("$D_3d$")
    plt.title("Fixed $A$ = %s" % Aval)
    figname = "D3d_fixed_A%i.png" % Aval
    plt.savefig(figname, bbox_inches="tight")
    print("Plot saved in %s." % figname)


args = docopt(__doc__)
outfile = args["<outfile>"]
print("Reading %s..." % outfile)
try:
    df = pd.read_csv(outfile, sep=" ", header=None)
except FileNotFoundError:
    sys.exit("File %s does not exist." % outfile)
df.columns = ["sys", "D3d"]
df["rd"] = [float(line.split("_")[1]) for line in df.sys]
df["A"] = [float(line.split("_")[3]) for line in df.sys]
df["B"] = [float(line.split("_")[5]) for line in df.sys]
rd = df.rd[0]
df = df[df.A < 0]              # process only negative As
#df.D3d[df.D3d < 0.0] = 0.0     # remove negative numbers
Dcut = float(args["--dcut"])   # lowest dr^2 = 0.015: D = 0.015/6 = 0.0025

df_lg = df[df.D3d >= Dcut].sys.values
outname = "configs_lg.out"
np.savetxt(outname, df_lg, fmt="%s")
print("Liquid and gas configs for cutoff %.3e saved in %s." % (Dcut, outname))

df_s = df[df.D3d < Dcut].sys.values
outname = "configs_s.out"
np.savetxt(outname, df_s, fmt="%s")
print("Solid configs for cutoff %.3e saved in %s." % (Dcut, outname))

df = df.drop("sys", 1)
df = df[["rd", "A", "B", "D3d"]]
dfname = "diff_cleared.out"
df.to_csv(dfname)
print("Cleared data frame saved in %s." % dfname)
df2 = df[["A", "B", "D3d"]]
df2.to_csv("diff_gnuplot.out", sep="\t", header=False)

As = np.array(sorted(list(set(df.A))))
Bs = np.array(sorted(list(set(df.B))))

if args["--plot_B_cut"]:
    for B in Bs:
        cut_df = df[["A", "D3d"]][df.B == B]
        cut_df = cut_df.sort_values(["A"])
        plot_B_cut(cut_df, B, rd)

if args["--plot_A_cut"]:
    for A in As:
        cut_df = df[["B", "D3d"]][df.A == A]
        cut_df = cut_df.sort_values(["B"])
        plot_A_cut(cut_df, A, rd)

if args["--plot-imshow"]:
    interp = args["--interp"]
    M, N = len(As), len(Bs)
    As_d = {k: v for k, v in zip(sorted(As), range(M))}
    Bs_d = {k: v for k, v in zip(sorted(Bs), range(N))}
    color_cut = 0.2
    df.D3d[df.D3d > color_cut] = color_cut        # set the colour scale
#    df.D3d[df.D3d > top_D3d[rd]] = top_D3d[rd]   # set the colour scale
    extent = np.array([5, 100, -5, -100]) + np.array([-1, 1, 1, -1]) * 2.5

    D3d = np.zeros((M, N))
    D3d_solid = np.zeros((M, N))
    for i, r in df.iterrows():
        if r["B"] > - r["A"] * 2 * np.pi * rd**3 / 15:
            D3d[As_d[r["A"]], Bs_d[r["B"]]] = r["D3d"]
            if r["D3d"] > Dcut:
                D3d_solid[As_d[r["A"]], Bs_d[r["B"]]] = 1.0
            else:
                D3d_solid[As_d[r["A"]], Bs_d[r["B"]]] = 0.0
        else:       # no-go region should be in white
            D3d[As_d[r["A"]], Bs_d[r["B"]]] = np.nan
            D3d_solid[As_d[r["A"]], Bs_d[r["B"]]] = np.nan
    plt.imshow(D3d, extent=extent, vmin=0.0, interpolation=interp)
    plt.xlabel("$B$")
    plt.ylabel("$A$")
    plt.colorbar()
    plt.savefig("diff_imshow.png", bbox_inches="tight")
    print("Heat map saved in diff_imshow.png.")

    plt.clf()
    print("Cutoff for solid D: %.5f" % Dcut)
    plt.imshow(D3d_solid, vmax=1.0, vmin=0.0, extent=extent, \
            interpolation=interp)
    plt.xlabel("$B$")
    plt.ylabel("$A$")
    plt.title("$D < %.5f$" % Dcut)
    plt.savefig("diff_imshow_solid.png", bbox_inches="tight")
    print("Heat map of solidified liquid saved in diff_imshow_solid.png.")


