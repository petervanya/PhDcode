#!/usr/bin/env python
"""Usage:
    ptfe_profile.py <files> (1d <axis> | 2d <plane> --depth <d> --thick <h>)
                     [--bins <bins> --subst <s> --boxsize <L>]

Create 1D or 2D profile of water molecules in Nafion w or wo electrodes.
* 1D: choose axis along which to profile
* 2D: choose plane and depth at which to take a slab of thickness <h>
      in DPD units (1 = 8.14 AA)

Arguments:
    <files>             Dump files from LAMMPS
    1d                  1d profile of water
    2d                  2d profile of water
    <axis>              Select the axis to profile along ("x", "y", "z")
    <plane>             Select the plane to profile on ("xy", "yz", "xz")

Options:
    --subst <s>         Requested substance: "water", "sulfonic", "backbone" [default: water]
    --depth <d>         Depth at which to profile the plane in DPD units [default: 20]
    --thick <h>         Thickness of the profile in DPD units [default: 1]
    --boxsize <L>       Box size [default: 40]
    --bins <bins>       Number of bins for the histogram [default: 150]

pv278@cam.ac.uk, 03/01/15
"""
import numpy as np
from math import log, sqrt
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import glob
import sys
from docopt import docopt
import lmp_lib as ll

rc = 8.14e-10

def read_outfile(outfile):
    """Read one xyz outfile into a numpy matrix"""
    A = open(outfile, "r").readlines()[2:]
    A = np.array([line.split() for line in A], order="F").astype(float)
    return A


def set_num_bins(N, method="sqrt"):
    """Set number of histogram bins, various recipes.
    Available methods: rice, sturges, sqrt (default)"""
    if method == "rice":
        return int(2*N**(1./3)) + 1   # Rice rule
    elif method == "sturges":
        return int(log(N, 2)+1) + 1   # Sturges' formula
    else:
        return int(sqrt(N)) + 1       # most primitive


def create_1d_profile(dumpfiles, axis, subst, bins):
    """NEW function with pre-set bins and fixed water plot"""
    Nfiles = len(dumpfiles)
    res = np.zeros(len(bins)-1)
    if subst == "water":
        for dumpfile in dumpfiles:
            A = read_outfile(dumpfile)
            beadsC = A[A[:, 0] == 3][:, 1:]   # bead C, 3 H2O molecules
            beadsW = A[A[:, 0] == 4][:, 1:]   # bead W, 6 H2O molecules
            profileC, bins = np.histogram(beadsC[:, axis], bins=bins)
            profileW, bins = np.histogram(beadsW[:, axis], bins=bins)
#            res += 3*profileC/float(Nfiles)
            res += 6*profileW/float(Nfiles)
    elif subst == "sulfonic":
        for dumpfile in dumpfiles:
            A = read_outfile(dumpfile)
            beads = A[A[:, 0] == 3][:, 1:]
            profile, bins = np.histogram(beads[:, axis], bins=bins)
            res += profile/float(Nfiles)
    elif subst == "backbone":
        for dumpfile in dumpfiles:
            A = read_outfile(dumpfile)
            beads1 = A[A[:, 0] == 1][:, 1:]   # bead A, 6 CF2 groups
            beads2 = A[A[:, 0] == 2][:, 1:]   # bead B, 5 CF2 groups
            beads3 = A[A[:, 0] == 3][:, 1:]   # bead B, 1 CF3 group
            profile1, bins = np.histogram(beads1[:, axis], bins=bins)
            profile2, bins = np.histogram(beads2[:, axis], bins=bins)
            profile3, bins = np.histogram(beads3[:, axis], bins=bins)
            res += 6*profile1/float(Nfiles)
            res += 5*profile2/float(Nfiles)
            res += profile3/float(Nfiles)
    return res, bins


def create_2d_profile(dumpfiles, plane, nbins, D, H):
    """2D histogram from given xyz files 
    at certain depth D and slab height H.
    SO FAR ONLY FOR WATER BEADS.
    """
    Nf = len(dumpfiles)
    res = np.zeros((nbins, nbins))
    axis = set([0, 1, 2]).difference(plane)
    axis = list(axis)[0]    # choose perpendicular axis to the given plane
    print("Perp. axis:", axis)
    for dumpfile in dumpfiles:
        A = read_outfile(dumpfile)
        A = A[(A[:, axis+1] > D-H/2) & (A[:, axis+1] < D+H/2)] # +1 to account for bead names col
        beads3 = A[A[:, 0] == 3][:, 1:]  # bead C, 3 H2O molecules
        beads4 = A[A[:, 0] == 4][:, 1:]  # bead W, 6 H2O molecules
        profile3, _, _ = np.histogram2d(beads3[:, plane[0]], beads3[:, plane[1]], bins=nbins)
        profile4, _, _ = np.histogram2d(beads4[:, plane[0]], beads4[:, plane[1]], bins=nbins)
        res += 3*profile3/float(Nf)
        res += 6*profile4/float(Nf)
    return res


if __name__ == "__main__":
    args = docopt(__doc__)
    L = float(args["--boxsize"])#*rc
    axes = {"x": 0, "y": 1, "z": 2}
    planes = {"xy": (0, 1), "yz": (1, 2), "xz": (0, 2)}
    subst_map = {"water": "W", "sulfonic": "S", "backbone": "B"}
    dumpfiles = glob.glob(args["<files>"])
    nbins = int(args["--bins"])
    subst = args["--subst"]
    if subst not in subst_map.keys():
        print("Choose substance from 'water', 'sulfonic', 'backbone'.")
        sys.exit()
    
    if not dumpfiles:
        print("No files captured, aborting.")
        sys.exit()

    print("===== Water profile =====")
    print("Substance:", subst, "| Bins:", nbins, "| Box size:", L, \
          "| xyz files:", len(dumpfiles))
    
    if args["1d"]:
        try:
            axis = axes[args["<axis>"]]
        except KeyError:
            print("Choose axis from 'x', 'y', 'z'.")
            sys.exit()
        bins = np.linspace(0, L, nbins+1)
        profile, bins = create_1d_profile(dumpfiles, axis, subst, bins)
        bins = bins[:-1] + np.diff(bins)/2.0

        outname = "profile_1d_%s_%if.out" % (subst_map[subst], len(dumpfiles))
        np.savetxt(outname, np.vstack((bins, profile)).T)
        print("Array saved in", outname)

    elif args["2d"]:
        try:
            plane = planes[args["<plane>"]]
        except KeyError:
            print("Choose plane from 'xy', 'yz', 'xz'.")
            sys.exit()
        D = float(args["--depth"])
        H = float(args["--thick"])
        print("2D profile | Plane: %s | Depth: %.1f | Thickness: %.1f") % \
             (args["<plane>"], D, H)

        profile = create_2d_profile(dumpfiles, plane, nbins, D, H)
        outname = "profile_2d.out"
        np.savetxt(outname, profile)
        print("Matrix saved in", outname)

        plt.imshow(profile, cmap="spectral")# plt.get_cmap("gray")) #cmap = cm.Greys_r)
        plt.axis("off")
        plotname = "profile_2d_" + args["<plane>"] + "_" + str(len(dumpfiles)) + "f.png"
        plt.savefig(plotname)
        print("Plot saved in", plotname)



