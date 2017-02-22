#!/usr/bin/env python
"""Usage:
    process_perc_bulk.py parse [--path <p>]
    process_perc_bulk.py plot <dfpath> [--sparse --pdf]

[AD HOC] process percolations of bulk Nafion. Collect all data about
percolation of water in Nafion into one table and plot relevant figures.

Arguments:
    parse       Create a master dataframe and save it
    <dfpath>    Path of the dataframe with diffusivities

Options:
    --path <p>     Path with logfiles [default: ./]
    --sparse       Process only water uptakes 6, 9, 12, 16
    --pdf          Export to pdf

pv278@cam.ac.uk, 12/02/17
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
    cols = ["lmbda", "P2d", "P3d"]
    df = DataFrame(index=range(len(lmbdas)), columns=cols)

    cnt = 0
    for l in lmbdas:
        data = [l]
        for dim in [2, 3]:
            fname = "clustering_%id_l%i.log" % (dim, l)
            fpath = str(default_path) + fname
            try:
                f = open(fpath, "r").readlines()
                data.append(float(f[-1].split()[-1]))
            except FileNotFoundError:
                print("File not found: %s." % fpath)

        if len(data) == 3:
            df.loc[cnt] = data
            cnt += 1
    return df


def plot(df, lmbdas, dim, to_pdf, is_sparse):
    """Plot 2d or 3d percolation
    * el: electrode, 'Carbon' or 'Quartz'
    * dim: '2d' or '3d'
    * lmbdas: range of lmbdas
    * is_sparse: plot only certain lmbdas
    """
    fig, ax = plt.subplots()    # plt.figure()
    Pinf = "P" + dim

    sel = np.array(df[["lmbda", Pinf]])
    print(sel)
    plt.plot(sel[:, 0], sel[:, 1], "D-", lw=4, ms=10, mew=0, label="Bulk")

    plt.xticks([4, 8, 12, 16, 20, 24])
    plt.xlim([lmbdas[0]-1, lmbdas[-1]+1])
    plt.xlabel("$\lambda$")
    plt.ylabel("Percolation cluster strength")
    plt.ylim([0, 1])
    plt.legend(loc=2, fontsize=22)

    fmt = "pdf" if to_pdf else "png"
    sp = "_sp" if is_sparse else ""
    plotname = "%s/Pinf_%s%s.%s" % (default_path, dim, sp, fmt)
    plt.savefig(plotname, bbox_inches='tight')


if __name__ == "__main__":
    args = docopt(__doc__)
    default_path = args["--path"]
    print("Directory: %s" % default_path)

#    lmbdas = np.arange(4, 25, 2).astype(int)
    lmbdas = sorted(list(range(4, 25, 2)) + [9])

    if args["parse"]:
        df = collect_data(default_path, lmbdas)
        dfname = default_path + "all_perc.csv"
        df.to_csv(dfname, delim=",")
        print("Dataframe saved in %s." % dfname)


    if args["plot"]:
        dfpath = args["<dfpath>"]
        if args["--sparse"]:
            lmbdas = [6, 9, 12, 16]
        else:
            lmbdas = sorted(list(range(4, 25, 2)))
        df = read_csv(dfpath)
        df = df.loc[df.lmbda.isin(lmbdas)]

        to_pdf = args["--pdf"]
        if to_pdf:
            plt.rc('text', usetex=True)
            plt.rc('font', family='serif')

        print("Plotting different electrodes.")
        print("2d percolation")
        plot(df, lmbdas, "2d", to_pdf, args["--sparse"])
        print("3d percolation")
        plot(df, lmbdas, "3d", to_pdf, args["--sparse"])


