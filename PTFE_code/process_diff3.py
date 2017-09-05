#!/usr/bin/env python
"""Usage: 
    process_diff3.py [--path <p> --pdf]

[AD HOC] Process diffusivities. Collect all the data about 
diffusivities of water in Nafion into one table 
and plot the relevant figures.

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
                        if "1d diff" in line:
                            data.extend(np.array(line.split()[5:]).astype(float))
                        if "2d diff" in line:
                            data.extend(np.array(line.split()[5:]).astype(float))
                        if "3d diff" in line:
                            data.extend(np.array(line.split()[2:]).astype(float))
                except FileNotFoundError:
                    print("File not found: %s." % fpath)

                if len(data) == len(cols):
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
    fig, ax = plt.subplots(figsize=(8, 6))

    for d in widths:
        lbl = "%i nm" % d
        sel = np.array(df[["lmbda", D]][(df.d == d) & (df.elect == el)])
        plt.plot(sel[:, 0], sel[:, 1], "+-", lw=2, ms=2, mew=2, label=lbl)

    plt.xlim([lmbdas[0]-1, lmbdas[-1]+1])
    plt.ylim([0, 0.25])
    plt.xticks([4, 8, 12, 16, 20, 24])
    plt.xlabel("$\lambda$")
    ylbl = {"Dx": "Normal", "Dyz": "Parallel", "D3d": "3d"}
#    plt.ylabel("$D_{\mathrm{%s}} \; (10^{-9} \; \mathrm{m^2/s}) $" % ylbl[D])
    plt.ylabel("$D_{\mathrm{%s}}$" % ylbl[D])

    plt.legend(loc=2, fontsize=22)
    plt.grid()
    plt.title(el)
    fmt = "pdf" if to_pdf else "png"
    plotname = "%s/%s_%s.%s" % (default_path, D, el.lower(), fmt)
    plt.savefig(plotname, bbox_inches='tight')


if __name__ == "__main__":
    args = docopt(__doc__)
    default_path = args["--path"]
    print("Directory: %s" % default_path)

    elects_d = {"c": "Carbon", "q": "Quartz"}
    elects = sorted(elects_d.values())

    lmbdas = list(range(4, 25, 2))
    widths = [5, 10, 15, 20]
    df = collect_data(default_path, elects_d, widths, lmbdas)
    dfname = default_path + "all_diff.csv"
    df.to_csv(dfname)
    print("Dataframe saved in %s." % dfname)

    to_pdf = args["--pdf"]
    if to_pdf:
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')

    print("Plotting diffusivities.")
    for el in elects:
        for diff in ["Dx", "Dyz", "D3d"]:
            print(el, diff)
            plot_el(df, el, diff, widths, lmbdas, to_pdf)


