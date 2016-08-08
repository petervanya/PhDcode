#!/usr/bin/env python
"""Usage:
    bandgap_vs_cluster.py (--spin | --gs) [-f <f>]

Plot bandgap vs number of atoms for lowest spin state

Options:
  --spin               Take spin = 0 state for each cluster
  --gs                 Take ground state for each cluster
  -f <f>,--format <f>  File format [default: png]

23/05/15
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from docopt import docopt


def get_Eg(cluster_list):
    Eg = []
    for cluster in cluster_list:
        inpath = base_dir + "/Pt" + cluster + "_bandgap.out"
        with open(inpath) as f:
            a = f.read().split()
        Eg.append(float(a[1]))
    return Eg


def get_cluster_list():
    """Read cluster configurations from the file"""
    with open("cluster_list.txt") as f:
        cluster_list = f.readlines()
    cluster_list = [i.strip("\n") for i in cluster_list]
    return cluster_list


def plot(Eg, atom_numbers, outpath, title):
    #plt.rc("text", usetex=True)
    #plt.rc("font", family="serif")
    fig = plt.figure()
    plt.plot(atom_numbers, Eg, "bo", ms=10)
    plt.xlabel("Number of atoms in cluster")
    plt.ylabel("$E_g$ (eV)")
    plt.title(title)
    
    fig.set_size_inches(4, 3)
    #plt.tight_layout()
    plt.savefig(outpath)
    print("Figure saved in",outpath)


if __name__ == "__main__":
    args = docopt(__doc__)
    cluster_list = get_cluster_list()
    atom_numbers = [sum([int(j) for j in i.split("_")]) for i in cluster_list]
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    base_dir += "/Plain/Outfiles/Evals/"
    fileformat = args["--format"]

    if args["--spin"]:
        Eg = get_Eg(cluster_list)
        outpath = base_dir + "bg_vs_cluster_S0." + fileformat
        plot(Eg, atom_numbers, outpath, "$S$ = 0")
    if args["--gs"]:
        A = [line.split() for line in open(base_dir + "bg_vs_cluster_gs.out").readlines()]
        A = np.array(A)
        print(A[:,2])
        Eg = [float(i) for i in A[:,2]]

        outpath = base_dir + "bg_vs_cluster_gs." + fileformat
        plot(Eg, atom_numbers, outpath, "Lowest spin state")






