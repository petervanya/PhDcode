#!/usr/bin/env python
"""Usage:
    process_coord_num.py <outfile> [--interp <in> --tex --fmt <fmt>]

Load the raw coord num file, create a data frame out of it and plot results.

Options:
    --interp <in>    Imshow interpolation [default: none]
    --tex            Use tex for axes
    --fmt <fmt>      Plot format [default: png]

13/11/17
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
df.columns = ["sys", "cn"]
df["rd"] = [float(line.split("_")[1]) for line in df.sys]
df["A"] = [float(line.split("_")[3]) for line in df.sys]
df["B"] = [float(line.split("_")[5]) for line in df.sys]
rd = df.rd[0]
df = df[df.A < 0]
df = df.drop("sys", 1)
df = df[["rd", "A", "B", "cn"]]

#As = np.array(sorted(list(set(df.A))))
#Bs = np.array(sorted(list(set(df.B))))
As = np.linspace(-5, -100, 20)
Bs = np.linspace(5, 100, 20)

interp = args["--interp"]
M, N = len(As), len(Bs)
As_d = {k: v for k, v in zip(sorted(As), range(M))}
Bs_d = {k: v for k, v in zip(sorted(Bs), range(N))}
extent = np.array([2.5, 102.5, -2.5, -102.5])

figfmt = args["--fmt"]
if args["--tex"]:
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

fig, ax = plt.subplots(figsize=(8, 6))
CN = np.empty((M, N))
CN[:] = np.nan   # initialise as nan to make white as default
nogo_A = np.arange(0, -110, -5)
nogo_B = -nogo_A * 2 * np.pi * rd**3 / 15

for _, r in df.iterrows():
    if r["B"] > - r["A"] * 2 * np.pi * rd**3 / 15:
#    if r["B"] > 25:
        CN[As_d[r["A"]], Bs_d[r["B"]]] = r["cn"]
    else:       # no-go region should be in white
        CN[As_d[r["A"]], Bs_d[r["B"]]] = np.nan

#if rd == 0.65:
#    CN[As_d[-100], Bs_d[15]] = np.nan
#if rd == 0.55:
#    CN[As_d[-100], Bs_d[10]] = np.nan

#plt.plot(nogo_B, nogo_A, "r", lw=3)
#xt_bcc = {0.55: 70, 0.65: 85, 0.75: 80}
#yt_bcc = {0.55: -60, 0.65: -70, 0.75: -90}
#xt_hex = {0.55: 40, 0.65: 35, 0.75: 30}
#yt_hex = {0.55: -90, 0.65: -95, 0.75: -95}

font_bcc = {"weight": "bold", "color": "white"}
font_hex = {"weight": "heavy", "color": "black"}
plt.imshow(CN, extent=extent, vmin=0, vmax=24, interpolation=interp)
#plt.text(xt_bcc[rd], yt_bcc[rd], "bcc", fontdict=font_bcc)
#plt.text(xt_hex[rd], yt_hex[rd], "hex", fontdict=font_hex)
plt.colorbar()
plt.xlabel("$B$")
plt.ylabel("$A$")
plt.title("$r_{\\rm{d}} = %s$" % rd)
figname = "cn_imshow_rd_%i_full." % (rd * 100) + figfmt
plt.savefig(figname, bbox_inches="tight")
print("Coord num heat map saved in %s." % figname)


