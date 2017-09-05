#!/usr/bin/env python
"""Usage:
    process_perc_bulk2.py [--path <p> --pdf]

[AD HOC] process percolations of bulk Nafion. Collect all data about
percolation of water in Nafion into one table and plot relevant figures.

Options:
    --path <p>     Path with logfiles [default: ./]
    --pdf          Export to pdf

pv278@cam.ac.uk, 04/05/17
"""
import numpy as np
from pandas import DataFrame, read_csv
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import matplotlib
from docopt import docopt

matplotlib.rcParams.update({'font.size': 24})


def collect_data(default_path, lmbdas):
    """Visit directory for each lmbda and
    parse the outfile for 2d and 3d percolations"""
    cols = ["lmbda", "P2d", "occ_2d", "P3d", "occ_3d"]
    df = DataFrame(index=range(len(lmbdas)), columns=cols)

    cnt = 0
    for l in lmbdas:
        data = [l]
        for dim in [2, 3]:
            fname = "clustering_%id_l%i.log" % (dim, l)
            fpath = str(default_path) + fname
            try:
                f = open(fpath, "r").readlines()
                data.append(float(f[-1].split()[-1]))  # P inf
                fields = [line for line in f if "Full fields" in line]
                data.append(fields[-1].split()[-1])
            except FileNotFoundError:
                print("File not found: %s." % fpath)

        if len(data) == len(cols):
            df.loc[cnt] = data
            cnt += 1
    return df


def plot(df, lmbdas, dim, to_pdf):
    """Plot 2d or 3d percolation
    * el: electrode, 'Carbon' or 'Quartz'
    * dim: '2d' or '3d'
    * lmbdas: range of lmbdas
    """
    fig, ax = plt.subplots(figsize=(8, 6))
    Pinf = "P" + dim

    sel = np.array(df[["lmbda", Pinf]])
    plt.plot(sel[:, 0], sel[:, 1], "D-", lw=4, ms=10, mew=0, label="Bulk")

    plt.xticks([4, 8, 12, 16, 20, 24])
    plt.xlim([lmbdas[0]-1, lmbdas[-1]+1])
    plt.xlabel("$\lambda$")
    plt.ylabel("$P_{\infty}$")
    plt.ylim([0, 1])
    plt.legend(loc=2, fontsize=22)

    fmt = "pdf" if to_pdf else "png"
    plotname = "%s/Pinf_%s.%s" % (default_path, dim, fmt)
    plt.savefig(plotname, bbox_inches='tight')


if __name__ == "__main__":
    args = docopt(__doc__)
    default_path = args["--path"]
    print("Processing percolation")
    print("Directory: %s" % default_path)

    lmbdas = np.arange(4, 25, 2).astype(int)
    df = collect_data(default_path, lmbdas)
    dfname = default_path + "all_perc.csv"
    df.to_csv(dfname)
    print("Dataframe saved in %s." % dfname)

    if df.empty:
        sys.exit("Error: Empty dataframe.")

    to_pdf = args["--pdf"]
    if to_pdf:
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')

    print("Plotting different electrodes.")
    print("2d percolation")
    plot(df, lmbdas, "2d", to_pdf)
    print("3d percolation")
    plot(df, lmbdas, "3d", to_pdf)


