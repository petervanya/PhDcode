#!/usr/bin/env python
"""Usage:
    process_diff_bulk.py --bt <bt> [--path <p> --sparse --pdf]

[AD HOC] Process diffusivities for bulk Nafion. Collect all data about
diffusivities of water in Nafion into one table 
and plot the relevant figures.

Arguments:
    --bt <bt>      Bead type

Options:
    --path <p>     Path with logfiles [default: ./]
    --sparse       Process only water uptakes 6, 9, 12, 16
    --pdf          Export to pdf

pv278@cam.ac.uk, 13/02/16
"""
import numpy as np
from pandas import DataFrame, read_csv
import matplotlib.pyplot as plt
import matplotlib
import sys
from docopt import docopt


matplotlib.rcParams.update({'font.size': 24})
bt_colors = {1: "red", 2: "red", 3: "green", 4: "blue", 5: "grey", 6: "brown"}


def collect_data(default_path, lmbdas, bt):
    """Visit directory for each lmbda and
    parse the outfile for 1d, 2d, and 3d diffusivities"""
    cols = ["lmbda", "Dx", "Dy", "Dz", "Dxy", "Dyz", "Dxz", "D3d"]
    df = DataFrame(index=range(len(lmbdas)), columns=cols)

    cnt = 0
    for l in lmbdas:
        data = [l]
        fname = "diffusivity_l%i_b%i.log" % (l, bt)
        fpath = default_path + fname
        try:
            f = open(fpath, "r").readlines()
            for line in f:
                if "1d diff" in line:
                    data.extend(np.array(line.split()[5:]).astype(float))
                if "2d diff" in line:
                    data.extend(np.array(line.split()[5:]).astype(float))
                if "3d diff" in line:
                    data.extend(np.array(line.split()[2:]).astype(float))
        except FileNotFoundError:
            print("File not found: %s." % fpath)

        if len(data) == 8:
            df.loc[cnt] = data
            cnt += 1
    return df


def plot(df, D, lmbdas, to_pdf, bt, is_sparse):
    """Plot parallel and normal diffusivities for electrode
    * df: dataframe
    * D: either 'Dx' or 'Dyz'
    * lmbdas: range of lmbdas
    * is_sparse: flag for filename
    """
    fig, ax = plt.subplots(figsize=(8, 6))
    ymax = 0.0

    sel = np.array(df[["lmbda", D]])
    plt.plot(sel[:, 0], sel[:, 1] * 1e9, \
            "D-", lw=4, ms=10, mew=0, label="Bulk", color=bt_colors[bt])

    ymax_temp = (max(sel[:, 1] * 1e9) // 10 + 1) * 10  # round to 10s
    if ymax_temp > ymax: ymax = ymax_temp
    plt.ylim([0.0, ymax])

    plt.xticks([4, 8, 12, 16, 20, 24])
    plt.xlim([lmbdas[0]-1, lmbdas[-1]+1])
    plt.yticks(np.arange(0.0, ymax+1, 5))
    plt.xlabel("$\lambda$")
    if D in ["Dx", "Dy", "Dz"]:
        ylbl = "Normal"
    elif D in ["Dxy", "Dyz", "Dxz"]:
        ylbl = "Parralel"
    elif D == "D3d":
        ylbl = "3d"
    plt.ylabel("$D_{\mathrm{%s}} \; (10^{-9} \; \mathrm{m^2/s}) $" % ylbl)

    plt.legend(loc="best", fontsize=22)

    fmt = "pdf" if to_pdf else "png"
    sp = "_sp" if is_sparse else ""
    plotname = "%s/%s%s_b%i.%s" % \
            (default_path, D, sp, bt, fmt)
    plt.savefig(plotname, bbox_inches='tight')


if __name__ == "__main__":
    args = docopt(__doc__)
    default_path = args["--path"]
    bt = int(args["--bt"])
    print("Processing diffusivty, bead type: %i" % bt)
    print("Directory: %s" % default_path)

    lmbdas = sorted(list(range(4, 25, 2)) + [9])
    df = collect_data(default_path, lmbdas, bt)
    dfname = default_path + "all_diff_b%i.csv" % bt
    df.to_csv(dfname)
    print("Dataframe saved in %s." % dfname)

    if args["--sparse"]:
        lmbdas = [6, 9, 12, 16]
    else:
        lmbdas = list(range(4, 25, 2))
    df = df.loc[df.lmbda.isin(lmbdas)]
    if df.empty:
        sys.exit("Error: Empty dataframe.")

    to_pdf = args["--pdf"]
    if to_pdf:           # use latex fonts
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')

    print("Plotting different electrodes.")
    print("Normal diffusivities.")
    plot(df, "Dx", lmbdas, to_pdf, bt, args["--sparse"])
    print("Parallel diffusivities.")
    plot(df, "Dyz", lmbdas, to_pdf, bt, args["--sparse"])
    print("Bulk diffusivities.")
    plot(df, "D3d", lmbdas, to_pdf, bt, args["--sparse"])


