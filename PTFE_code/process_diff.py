#!/usr/bin/env python
"""Usage: 
    process_diff.py parse
    process_diff.py plot by_elect <dfpath> [--pdf]
    process_diff.py plot by_lambda <dfpath> [--pdf]

[AD HOC] Process diffusivities. Collect all the data about 
diffusivities of water in Nafion into one table 
and plot the relevant figures.

Arguments:
    parse       Create a master dataframe and save it
    by_elect    Create plots of Dx and Dyz for different electrodes
    by_lambda   Create plots of Dx and Dyz for different water uptakes
    <dfpath>    Path of the dataframe with diffusivities

Options:
    --pdf       Export in pdf, does not work on Cottrell

pv278@cam.ac.uk, 17/11/16
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import matplotlib
import pandas as pd
from pandas import DataFrame
from docopt import docopt


matplotlib.rcParams.update({'font.size': 16})


def collect_data(default_path, elects_d, widths, lmbdas):
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
                f = open(fpath, "r").readlines()
                for line in f:
                    if "1d in SI" in line:
                        data.extend(np.array(line.split()[3:]).astype(float))
                    if "2d in SI" in line:
                        data.extend(np.array(line.split()[3:]).astype(float))
                    if "3d in SI" in line:
                        data.extend(np.array(line.split()[3:]).astype(float))
                df.loc[cnt] = data
                cnt += 1
    return df


def plot_el(df, el, D, widths, lmbdas, to_pdf):
    """Plot parallel and normal diffusivities for electrode
    * df: dataframe
    * el: electrode
    * D: either 'Dx' or 'Dyz'
    * widths: range of widths
    * lmbdas: range of lmbdas
    """
    df_cut = df[df.elect == el]
    fig, ax = plt.subplots()    # plt.figure()
#    plt.tight_layout()
    plt.gcf().subplots_adjust(left=0.18)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
    for l in lmbdas:
        sel = df_cut[["d", D]][df_cut.lmbda == l]
        sel = np.array(sel)
        lbl = "$\lambda = %i$" % l
        plt.plot(sel[:, 0], sel[:, 1], "+-", lw=2, ms=10, mew=2, label=lbl)
        plt.xticks(widths)
        plt.xlim([widths[0]-1, widths[-1]+1])
        plt.xlabel("Electrode width (nm)")
        ylbl = "Normal" if len(D) == 2 else "Parralel"
        plt.ylabel("$D_{\mathrm{%s}}$" % ylbl)
        plt.legend(loc="best")
        plt.title(el)
    fmt = ".pdf" if to_pdf else ".png"
    plotname = default_path + D + "_" + el.lower() + fmt
    plt.savefig(plotname)


def plot_lmbda(df, lmbda, D, widths, elects, to_pdf):
    """Plot parallel and normal diffusivities for a given uptake
    * df: dataframe
    * lmbda: water uptake
    * D: either 'Dx' or 'Dyz'
    * widths: range of widths
    * elects: list of electrodes
    """
    df_cut = df[df.lmbda == lmbda]
    fig, ax = plt.subplots()    # plt.figure()
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
    for el in elects:
        sel = df_cut[["d", D]][df_cut.elect == el]
        sel = np.array(sel)
        plt.plot(sel[:, 0], sel[:, 1], "+-", lw=2, ms=10, mew=2, label=el)
        plt.xticks(widths)
        plt.xlim([widths[0]-1, widths[-1]+1])
        plt.xlabel("Electrode width (nm)")
        ylbl = "Normal" if len(D) == 2 else "Parralel"
        plt.ylabel("$D_{\mathrm{%s}}$" % ylbl)
        plt.legend(loc="best")
        plt.title("$\lambda = %i$" % lmbda)
    fmt = ".pdf" if to_pdf else ".png"
    plotname = default_path + D + "_l" + str(lmbda) + fmt
    plt.savefig(plotname)


if __name__ == "__main__":
    args = docopt(__doc__)
    default_path = "/home/pv278/Res/PTFEsim/Slab_DLMS/DataCommon/"

    lmbdas = [6, 9, 12]
    widths = [5, 10, 15]
    elects_d = {"c": "Carbon", "q": "Quartz", "s": "Silica"}
    elects = sorted(elects_d.values())

    if args["parse"]:
        df = collect_data(default_path, elects_d, widths, lmbdas)
        dfname = default_path + "all_diff.csv"
        df.to_csv(dfname, delim=",")
        print("Dataframe saved in %s." % dfname)

    if args["plot"]:
        dfpath = args["<dfpath>"]
        df = pd.read_csv(dfpath)
        to_pdf = args["--pdf"]
        if to_pdf:   # use latex for everything
            plt.rc('text', usetex=True)
            plt.rc('font', family='serif')

        if args["by_elect"]:
            print("Plotting different electrodes.")
            print("Normal diffusivities.")
            for el in elects:
                print("Plotting %s" % el)
                plot_el(df, el, "Dx", widths, lmbdas, to_pdf)
 
            print("Parallel diffusivities.")
            for el in elects:
                print("Plotting %s" % el)
                plot_el(df, el, "Dyz", widths, lmbdas, to_pdf)
 
        if args["by_lambda"]:
            print("Plotting different water uptakes.")
            print("Normal diffusivities.")
            for l in lmbdas:
                print("Plotting l = %s" % l)
                plot_lmbda(df, l, "Dx", widths, elects, to_pdf)
 
            print("Parallel diffusivities.")
            for l in lmbdas:
                print("Plotting l = %s" % l)
                plot_lmbda(df, l, "Dyz", widths, elects, to_pdf)


