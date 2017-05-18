#!/usr/bin/env python
"""Usage:
    process_perc.py parse [--path <p>]
    process_perc.py plot <dfpath> [--path <p> --sparse --pdf]

[AD HOC] process percolations. Collect all data about
percolation of water in Nafion into one table and
plot relevant figures.

Arguments:
    parse          Create a master dataframe and save it
    <dfpath>       Path of the dataframe with diffusivities

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


def collect_data(default_path, elects_d, widths, lmbdas):
    """Visit directory for each slab width and lmbda and
    parse the outfile for 2d and 3d percolations"""
    cols = ["elect", "d", "lmbda", "P2d", "P3d"]
    index = np.arange(len(lmbdas) * len(widths) * len(elects_d))
    df = DataFrame(index=index, columns=cols)

    cnt = 0
    for el in sorted(elects_d.keys()):
        for d in widths:
            for l in lmbdas:
                data = [elects_d[el], d, l]
                for dim in [2, 3]:
                    fname = "clustering_%id_%s_d%i_l%i.log" % (dim, el, d, l)
                    fpath = str(default_path) + fname
                    try:
                        f = open(fpath, "r").readlines()
                        data.append(float(f[-1].split()[-1]))
                    except FileNotFoundError:
                        print("File not found: %s." % fpath)

                if len(data) == 5:
                    df.loc[cnt] = data
                    cnt += 1
    return df


def plot_el(df, el, widths, lmbdas, dim, to_pdf, is_sparse):
    """Plot 2d or 3d percolation
    * df: dataframe
    * el: electrode, 'Carbon' or 'Quartz'
    * dim: '2d' or '3d'
    * widths: range of widths
    * lmbdas: range of lmbdas
    * is_sparse: plot only certain lmbdas
    """
    df_cut = df[df.elect == el]
    fig, ax = plt.subplots()    # plt.figure()
    Pinf = "P" + dim

    for d in widths:
        lbl = "%i nm" % d
        sel = np.array(df_cut[["lmbda", Pinf]][df_cut.d == d])
        plt.plot(sel[:, 0], sel[:, 1], "D-", lw=4, ms=10, mew=0, label=lbl)
   
    plt.xlabel("$\lambda$")
    plt.ylabel("Percolation cluster strength")
    plt.xticks([4, 8, 12, 16, 20, 24])
    plt.xlim([lmbdas[0]-1, lmbdas[-1]+1])
    plt.ylim([0, 1])
    plt.legend(loc=2, fontsize=22)
    plt.title(el)
    fmt = "pdf" if to_pdf else "png"
    sp = "_sp" if is_sparse else ""
    plotname = "%s/Pinf_%s_%s%s.%s" % \
            (default_path, dim, el.lower(), sp, fmt)
    plt.savefig(plotname, bbox_inches='tight')


if __name__ == "__main__":
    args = docopt(__doc__)
    default_path = args["--path"]
    print("Directory: %s" % default_path)

    lmbdas = sorted(list(range(4, 25, 2)) + [9])
    widths = np.array([5, 10, 15, 20])
    elects_d = {"c": "Carbon", "q": "Quartz"}
    elects = sorted(elects_d.values())

    if args["parse"]:
        df = collect_data(default_path, elects_d, widths, lmbdas)
        dfname = default_path + "all_perc.csv"
        df.to_csv(dfname, delim=",")
        print("Dataframe saved in %s." % dfname)


    if args["plot"]:
        if args["--sparse"]:
            lmbdas = [6, 9, 12, 16]
            widths = [5, 10, 15]
        else:
            lmbdas = list(range(4, 25, 2))
            widths = [5, 10, 15, 20]
        dfpath = args["<dfpath>"]
        df = read_csv(dfpath)
        df = df.loc[df.lmbda.isin(lmbdas) & df.d.isin(widths)]

        to_pdf = args["--pdf"]
        if to_pdf:
            plt.rc('text', usetex=True)
            plt.rc('font', family='serif')

        print("Plotting different electrodes.")
        print("2d percolation")
        for el in elects:
            print("Plotting %s" % el)
            plot_el(df, el, widths, lmbdas, "2d", to_pdf, args["--sparse"])
        print("3d percolation")
        for el in elects:
            print("Plotting %s" % el)
            plot_el(df, el, widths, lmbdas, "3d", to_pdf, args["--sparse"])


