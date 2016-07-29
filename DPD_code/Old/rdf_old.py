#!/usr/bin/env python
"""Usage:
    rdf.py gen_data <fnames> [--beads <b>] [--binsize <bs>]
    rdf.py plot hist <fnames> [--picname <pn>]
    rdf.py plot deriv2 <fnames> [--smooth <sm>]

Read all LAMMPS data files in the directory and compute the pair correlation function.
Uses Fortran routines from mat_ops.f95 produced with f2py.
* gen_data  generate the histogram of distances and save it in rdf.out
* plot      plot rdf.out
WARNING: PBC ARE NOT CORRECTLY IMPLEMENTED

Arguments:
    <fnames>         Regex for all the required xyz files

Options:
    --binsize <bs>   Size of the bins for the histogram, CAREFUL WITH UNITS [default: 0.05]
    --beads <b>      Bead type, 'all' or 1,...,n [default: all]
    --smooth <sm>    Smooth the data with a given span [default: 5]

pv278@cam.ac.uk, 02/09/15
"""
import numpy as np
from numpy.linalg import norm
import scipy.ndimage as ndimage
import matplotlib.pyplot as plt
from math import pi
import os, glob, sys
from docopt import docopt
import mat_ops           # Fortran .so file

bins = np.arange(0)

def read_outfile(outfile):
    """Read one xyz outfile into a numpy matrix"""
    A = open(outfile, "r").readlines()[2:]
    A = [line.split() for line in A]
    A = np.array(A, order="F").astype(float)
    return A


def save_data(outfile, *args):
    """Save two vectors into file"""
    m, n = len(args), len(args[0])   # cols, rows
    args = zip(*args)
    with open(outfile, "w") as f:
        for i in range(n):
            line = ""
            for j in range(m):
                line += str(args[i][j]) + "\t"
            line += "\n"
            f.write(line)


def plot_hist(r, vals, picname="hist.png", title="title", which="all", ymax=1000):
    """Plot a rdf function from given values"""
    plt.clf()
    plt.plot(r, vals)
    plt.ylim([0, ymax])
    plt.xlabel("$r$ (DPD units)")
    plt.ylabel("$g(r)$")
    plt.title(title)
    plt.savefig(picname)
    print("rdf saved into", picname)


def plot_deriv2(r, vals, picname="deriv2.png", title="title", which="all"):
    """Plot a 2nd derivative rdf function"""
    plt.clf()
    plt.plot(r, vals)
    plt.xlabel("$r$ (DPD units)")
    plt.title(title)
    plt.savefig(picname)
    print("2nd derivative saved into", picname)


def smooth(data, span=5):
    """Same as MATLAB smooth(data)"""
    N = len(data)
    data2 = np.zeros(N)
    n = (span-1)/2
    for i in range(n):
        data2[i] = sum(data[0:2*i+1])/(2*i+1)
        data2[N-1-i] = sum(data[N-(2*i+1):N])/(2*i+1)
    for i in range(n, N-n):
        data2[i] = sum(data[i-n:i+n+1])/span
    return data2


def deriv2(data, h=0.1):
    """Return array of 2nd derivatives"""
    N = len(data)
    data2 = np.zeros(N)
    data2[0] = (data[2] - 2*data[1] + data[0])/h**2
    data2[N-1] = (data[N-1] - 2*data[N-2] + data[N-3])/h**2
    for i in range(1, N-1):
        data2[i] = (data[i+1] - 2*data[i] + data[i-1])/h**2
    return data2


def get_one_rdf(outfile, dr=0.05, bead_type="all"):
    """Compute radial dist'n fcn from the xyz matrix using Fortran routine"""
    global bins
    A = read_outfile(outfile)
    if bead_type == "all":
        xyz_mat = A[:, 1:]
    else:
        xyz_mat = A[A[:, 0] == int(bead_type)][:, 1:]
    max_num = np.max(xyz_mat)
    if dr > max_num:
        print("WARNING: Bin size larger that box size: %e vs %e" % (dr, max_num))
    N = len(xyz_mat)
    d = mat_ops.get_pair_dist(xyz_mat)        # call Fortran
    print("Fortran routine done.")
    if len(bins) == 0:
        bins = np.arange(0, max_num+dr, dr)       # take box size as max entry of A
    rdf_raw, r = np.histogram(d, bins=bins)
    r = r[:-1] + np.diff(r)/2.0
    rdf = rdf_raw/(4*pi*r**2*dr)
    return r, rdf


def get_master_rdf(outfiles, dr=0.05, bead_type="all"):
    """Construct an rdf from all the available xyz files"""
    Nfiles = len(outfiles)
    rdf_mat = []
    for outfile in outfiles:
        r, rdf_i = get_one_rdf(outfile, dr, bead_type)
        rdf_mat.append(rdf_i)
        print(outfile, "done.")
    rdf_mat = np.array(rdf_mat).T
    np.savetxt("rdf_mat.out", rdf_mat)
    print("rdf matrix saved in rdf_mat.out")
    rdf = np.array(np.sum(rdf_mat, 1)/len(outfiles))
    return r, rdf


if __name__ == "__main__":
    args = docopt(__doc__)
#    print args
    if args["gen_data"]:
        dr = float(args["--binsize"])
        bead_type = args["--beads"]
        outfiles = glob.glob(args["<fnames>"])
        if len(outfiles) == 0:
            print("ERROR: No xyz files found. Aborting.")
            sys.exit()
        print(outfiles)
        Nfiles = len(outfiles)
        N = int(open(outfiles[0], "r").readline())
        print("Total particles:", N)
        r, vals = get_master_rdf(outfiles, dr, bead_type)
        fname = "rdf_" + bead_type + ".out"
        save_data(fname, r, vals)
        print("rdf saved in", fname)

    elif args["plot"] and args["hist"]:   # reads rdf.out and creates histograms
        outfiles = glob.glob(args["<fnames>"])
        for outfile in outfiles:
            A = np.loadtxt(outfile)
            r, vals = A[:, 0], A[:, 1]
            outfile = outfile[:-4]
            bead_type = outfile.split("_")[1]    # "all", "1", etc
            fA = int(outfile.split("_")[3])
            if bead_type == "1":                  # set ymax in advance
                ymax = fA/10.0*4000
            elif bead_type == "2":    # OPTIMISE THIS!
                ymax = (10.0-fA)/10.0*4000
            elif bead_type == "all":
                ymax = 10000

            dirname = "_".join(outfile.split("_")[2:])
            picname = "hist_" + bead_type + "_" + dirname + ".png"
            title = "beads " + bead_type
#            span = int(args["--smooth"])
#            assert (span-1)%2 == 0     # has to be odd number
#            vals = smooth(vals, span=span)
            plot_hist(r, vals, picname, title=title, which=bead_type, ymax=ymax)

    elif args["plot"] and args["deriv2"]:   # reads rdf.out and plots 2nd derivative of histograms
        outfiles = glob.glob(args["<fnames>"])
        print(outfiles)
        for outfile in outfiles:
            A = np.loadtxt(outfile)
            r, vals = A[:, 0], A[:, 1]
            outfile = outfile[:-4]
            bead_type = outfile.split("_")[1]    # "all", "1", etc
            dirname = "_".join(outfile.split("_")[2:])
            picname = "deriv2_" + bead_type + "_" + dirname + ".png"
            title = "beads " + bead_type

            span = int(args["--smooth"])
            assert (span-1)%2 == 0     # has to be odd number
            vals = smooth(vals, span=span)

            h = r[1] - r[0]
            vals2der = deriv2(vals, h)

            N = len(r)
            plot_deriv2(r[N/20:], vals2der[N/20:], picname, title=title, which=bead_type)
        

