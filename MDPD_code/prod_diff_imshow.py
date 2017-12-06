#!/usr/bin/env python
"""Usage:
    prod_diff_imshow.py <outfile> [--dcut <dc> --tex --pdf]
                        [--plot-imshow --interp <in>]
                        [--tex --fmt <fmt>]

Load the raw density file, create a data frame out of it and plot results.

Options:
    --plot-imshow    Plot 3D image
    --dcut <dc>      Cutoff for diffusivity of a solid [default: 0.0025]
    --interp <in>    Imshow interpolation [default: none]
    --tex            Use tex for axes
    --fmt <fmt>      Plot format [default: png]

08/09/17, modified 23/10/17
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import sys
from docopt import docopt

matplotlib.rcParams.update({'font.size': 20})


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
Dcut = float(args["--dcut"])   # lowest dr^2 = 0.015: D = 0.015/6 = 0.0025
df = df.drop("sys", 1)
df = df[["rd", "A", "B", "D3d"]]

As = np.array(sorted(list(set(df.A))))
Bs = np.array(sorted(list(set(df.B))))

interp = args["--interp"]
M, N = len(As), len(Bs)
As_d = {k: v for k, v in zip(sorted(As), range(M))}
Bs_d = {k: v for k, v in zip(sorted(Bs), range(N))}
color_cut = 0.2
df.D3d[df.D3d > color_cut] = color_cut        # set the colour scale
extent = np.array([2.5, 102.5, -2.5, -102.5])

fig, ax = plt.subplots(figsize=(8, 6))
figfmt = args["--fmt"]
if args["--tex"]:
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

D3d = np.zeros((M, N))
nogo_A = np.arange(0, -110, -5)
nogo_B = -nogo_A * 2 * np.pi * rd**3 / 15

for _, r in df.iterrows():
    D3d[As_d[r["A"]], Bs_d[r["B"]]] = r["D3d"]
plt.plot(nogo_B, nogo_A, "r", lw=3)
#ax.fill_between(nogo_B, nogo_A, where=nogo_B>=100, color="red")
plt.imshow(D3d, extent=extent, vmin=0.0, interpolation=interp)
plt.xlabel("$B$")
plt.ylabel("$A$")
plt.title("$r_{\\rm{d}} = %s$" % rd)
plt.colorbar()
figname = "diff_imshow." + figfmt
plt.savefig(figname, bbox_inches="tight")
print("Heat map saved in %s." % figname)

plt.clf()
print("Cutoff for solid D: %.5f" % Dcut)
plt.plot(nogo_B, nogo_A, "r", lw=3)
plt.imshow(D3d > Dcut, vmax=1.0, vmin=0.0, extent=extent, interpolation=interp)
plt.xlabel("$B$")
plt.ylabel("$A$")
plt.title("$r_{\\rm{d}} = %s$, $D < %.5f$" % (rd, Dcut))
figname = "diff_imshow_solid." + figfmt
plt.savefig(figname, bbox_inches="tight")
print("Heat map of solidified liquid saved in %s." % figname)


