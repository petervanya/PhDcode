#!/usr/bin/env python
"""
Analyse the surface tensions of various
spin configurations of Pt clusters w or w/o water

Usage:
    get_surf_tension.py <cluster> -e <e> [-d <d>] (--print | --save [--latex])
                        [--ext <ext>]

Arguments:
    <cluster>          Pt cluster, e.g. 9_10_9
Options:
    -d <d>,--dir <d>   Custom directory w input files (TODO)
    -e <e>,--eta <e>   Adsorption site, 1, 2 or 3
    --print            Print table on screen
    --save             Save table in ~/Platinum/Water/Outfiles
    --ext <ext>        File extension, e.g. "nosymm"
    --latex            Add latex formatting for output

pv278@cam.ac.uk, 07/04/15
"""
from docopt import docopt
import numpy as np
from math import sqrt
import pandas as pd
import os
import sys
from iolib import read_table, print_table, save_table


def isfloat(value):
    """check if character is float"""
    try:
        float(value)
        return True
    except ValueError:
        return False

    
def parse_data(file):
    data = []
    f = open(file,"r")
    for line in f:
        data.append( line.rstrip("\n").split(" \t ") )
    f.close()
    return data


if __name__ == "__main__":
    args = docopt(__doc__,version=1.0)
#    print args

    cluster = args["<cluster>"]
    eta = int(args["--eta"])
    ext = args["--ext"]

    fin = "Pt" + cluster + "_summary.out"
    Pt_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    inpath_p = Pt_dir + "/Plain/Outfiles/" + fin
    inpath_w = Pt_dir + "/Water/Outfiles/" + fin

    Ewater = -76.393551              # ZP corrected B3LYP/LANL2DZ
    a = 2.775
    S = a**2*sqrt(3.0)/2.0*1e-20     # in m**2
    Nspins = 11
    spin_list = range(Nspins)
    
#    names = ["spin", "converged", "E", "cycles", "error", "time"]
    Ap = pd.read_table(inpath_p, header=None, index_col=False)
    if ext:
        inpath_w += "." + ext
    Aw = pd.read_table(inpath_w, header=None)[Nspins*(eta-1) : Nspins*eta]
    Aw.index = spin_list

    res = pd.DataFrame(columns=["$E_\mathrm{b}$ (eV)", "$\sigma$ (J/m$^2$)"])
    for i in spin_list:
        Ep = Ap.ix[i][2]
        Ew = Aw.ix[i][6]
        if pd.notnull(Ew) and pd.notnull(Ep):
            dE = float(Ew) - Ewater - float(Ep)
            dE *= 27.21138505
            res.loc[spin_list[i]] = [dE, dE*1.602e-19/S]
    res.index.name = "Spin"
    
    if args["--print"]:
         print res
    if args["--save"]:
        fout = "Pt" + cluster + "_E" + str(eta) + "_sigma.out"
        outpath = Pt_dir + "/Water/Outfiles/" + fout
        if args["--ext"]:
            outpath += "." + ext
        if args["--latex"]:
            res.to_latex(outpath, escape=False) #outpath)
        else:
            res.to_csv(outpath, sep="\t")
        print "Table saved in",outpath
    
    
