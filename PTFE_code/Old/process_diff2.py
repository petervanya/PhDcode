#!/usr/bin/env python
"""Usage: 
    process_diff2.py parse [--path <p>]
    process_diff2.py plot <dfpath> [--path <p> --sparse --pdf]

[AD HOC] Process diffusivities. Collect all the data about 
diffusivities of water in Nafion into one table 
and plot the relevant figures.

Arguments:
    parse          Create a master dataframe and save it
    <dfpath>       Path of the dataframe with diffusivities

Options:
    --path <p>     Path with logfiles [default: ./]
    --sparse       Process only water uptakes 6, 9, 12, 16
    --pdf          Export to pdf

pv278@cam.ac.uk, 13/02/17
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
    parse the outfile for 1d, 2d, and 3d diffusivities"""
    cols = ["elect", "d", "lmbda", \
            "Dx", "Dy", "Dz", "Dxy", "Dyz", "Dxz", "D3d"]
    index = np.arange(len(lmbdas) * len(widths) * len(elects_d))
    df = DataFrame(index=index, columns=cols)

    cnt = 0
    for el in sorted(elects_d.keys()):
        for d in widths:
            for l in lmbdas:
                data = [elects_d[el], d, l]
                fname = "diffusivity_%s_d%i_l%i.log" % (el, d, l)
                fpath = default_path + fname
                try:
                    f = open(fpath, "r").readlines()
                    for line in f:
                        if "1d in SI" in line:
                            data.extend(np.array(line.split()[3:]).astype(float))
                        if "2d in SI" in line:
                            data.extend(np.array(line.split()[3:]).astype(float))
                        if "3d in SI" in line:
                            data.extend(np.array(line.split()[3:]).astype(float))
                except FileNotFoundError:
                    print("File not found: %s." % fpath)

                if len(data) == 10:
                    df.loc[cnt] = data
                    cnt += 1
    return df


def plot_el(df, el, D, widths, lmbdas, to_pdf, is_sparse):
    """Plot parallel and normal diffusivities for electrode
    * df: dataframe
    * el: electrode
    * D: either 'Dx' or 'Dyz'
    * widths: range of widths
    * lmbdas: range of lmbdas
    * is_sparse: flag for filename
    """
    df_cut = df[df.elect == el]
    fig, ax = plt.subplots(figsize=(8, 6))
    ymax = 0.0

    for d in widths:
        lbl = "%i nm" % d
        sel = np.array(df_cut[["lmbda", D]][df_cut.d == d])

        if el == "Carbon":
            Nl = len(lmbdas)
            eb = [1.0/i for i in range(1, Nl+1)]
            plt.errorbar(sel[:, 0], sel[:, 1] * 1e9, yerr=eb, \
                    fmt="D-", lw=4, ms=10, mew=0, label=lbl)
        else:
            plt.plot(sel[:, 0], sel[:, 1] * 1e9, \
                    "D-", lw=4, ms=10, mew=0, label=lbl)

        ymax_temp = (max(sel[:, 1] * 1e9) // 10 + 1) * 10  # round to 10s
        if ymax_temp > ymax: ymax = ymax_temp
        plt.ylim([0.0, ymax])

    plt.xlim([lmbdas[0]-1, lmbdas[-1]+1])
    plt.xticks([4, 8, 12, 16, 20, 24])
    plt.yticks(np.arange(0.0, ymax+1, 5))
    plt.xlabel("$\lambda$")
    ylbl = "Normal" if len(D) == 2 else "Parralel"
    plt.ylabel("$D_{\mathrm{%s}} \; (10^{-9} \; \mathrm{m^2/s}) $" % ylbl)

    plt.legend(loc="best", fontsize=22)
#    if el == "Carbon" and D == "Dx":
#        plt.legend(loc=2, fontsize=22)
#    if el == "Carbon" and D == "Dyz":
#        plt.legend(loc=(0.5, 0.5), fontsize=22)

    plt.grid()
    plt.title(el)
    fmt = "pdf" if to_pdf else "png"
    sp = "_sp" if is_sparse else ""
    plotname = "%s/%s_%s%s.%s" % \
            (default_path, D, el.lower(), sp, fmt)
#    plotname = default_path + D + "_" + el.lower() + fmt
    plt.savefig(plotname, bbox_inches='tight')


if __name__ == "__main__":
    args = docopt(__doc__)
    default_path = args["--path"]
    print("Directory: %s" % default_path)

    elects_d = {"c": "Carbon", "q": "Quartz"}
    elects = sorted(elects_d.values())

    if args["parse"]:
        lmbdas = sorted(list(range(4, 25, 2)) + [9])
        widths = np.array([5, 10, 15, 20])
        df = collect_data(default_path, elects_d, widths, lmbdas)
        dfname = default_path + "all_diff.csv"
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
        if to_pdf:           # use latex fonts
            plt.rc('text', usetex=True)
            plt.rc('font', family='serif')

        print("Plotting different electrodes.")
        print("Normal diffusivities.")
        for el in elects:
            print("Plotting %s" % el)
            plot_el(df, el, "Dx", widths, lmbdas, to_pdf, args["--sparse"])
        print("Parallel diffusivities.")
        for el in elects:
            print("Plotting %s" % el)
            plot_el(df, el, "Dyz", widths, lmbdas, to_pdf, args["--sparse"])


