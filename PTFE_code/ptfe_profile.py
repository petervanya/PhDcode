#!/usr/bin/env python
"""Usage:
    ptfe_profile.py 1d <ax> <files> [--subst <s> --L <L> --nbins <nb>]
    ptfe_profile.py 2d <pl> <files> [--depth <d> --thick <h> --L <L> --nbins <nb>]

Create 1D or 2D profile of water molecules in Nafion w or wo electrodes.
Beads: A, B, C (3 H2O), W (6 H2O)
* 1D: choose axis along which to profile
* 2D: choose plane and depth at which to take a slab of thickness <h>
      in DPD units (1 = 8.14 AA)

Arguments:
    <files>        xyz files
    <ax>           Axis to profile along ("x", "y", "z")
    <pl>           Plane to profile on ("xy", "yz", "xz")

Options:
    --subst <s>    Substance: "water", "sulfonic", "backbone", "el" [default: water]
    --depth <d>    Relative depth at which to profile the plane [default: 0.5]
    --thick <h>    Thickness of the profile in DPD units [default: 1]
    --L <L>        Box size [default: 40]
    --nbins <nb>   Number of bins for the histogram [default: 100]

pv278@cam.ac.uk, 03/01/15, modified 09/11/16
"""
import numpy as np
from math import sqrt
import matplotlib.pyplot as plt
import glob, sys, time
from docopt import docopt

rc = 8.14e-10

def read_outfile(outfile):
    """Read one xyz outfile into a numpy matrix"""
    A = open(outfile, "r").readlines()[2:]
    A = np.array([line.split() for line in A], order="F").astype(float)
    return A


def create_1d_profile(dumpfiles, axis, subst, bins):
    """Function with pre-set bins and fixed water plot"""
    Nfiles = len(dumpfiles)
    res = np.zeros(len(bins)-1)
    if subst == "water":
        for dumpfile in dumpfiles:
            A = read_outfile(dumpfile)
            beadsC = A[A[:, 0] == 3][:, 1:]   # bead C, 3 H2O molecules
            beadsW = A[A[:, 0] == 4][:, 1:]   # bead W, 6 H2O molecules
            profileC, bins = np.histogram(beadsC[:, axis], bins=bins)
            profileW, bins = np.histogram(beadsW[:, axis], bins=bins)
            res += 0.5 * profileC / float(Nfiles)
            res += 1.0 * profileW / float(Nfiles)
    elif subst == "sulfonic":
        for dumpfile in dumpfiles:
            A = read_outfile(dumpfile)
            beads = A[A[:, 0] == 3][:, 1:]
            profile, bins = np.histogram(beads[:, axis], bins=bins)
            res += 0.25 * profile / float(Nfiles)
    elif subst == "backbone":
        for dumpfile in dumpfiles:
            A = read_outfile(dumpfile)
            beads1 = A[A[:, 0] == 1][:, 1:]   # bead A, 6 CF2 groups
            beads2 = A[A[:, 0] == 2][:, 1:]   # bead B, 5 CF2 groups
            beads3 = A[A[:, 0] == 3][:, 1:]   # bead B, 1 CF3 group
            profile1, bins = np.histogram(beads1[:, axis], bins=bins)
            profile2, bins = np.histogram(beads2[:, axis], bins=bins)
            profile3, bins = np.histogram(beads3[:, axis], bins=bins)
            res += 1.0 * profile1 / float(Nfiles)
            res += 1.0 * profile2 / float(Nfiles)
            res += 0.25 * profile3 / float(Nfiles)
    elif subst == "el":
        for dumpfile in dumpfiles:
            A = read_outfile(dumpfile)
            beads = A[A[:, 0] == 5][:, 1:]
            profile, bins = np.histogram(beads[:, axis], bins=bins)
            res += profile / float(Nfiles)
    return res, bins


def create_2d_profile(dumpfiles, plane, nbins, D, H):
    """2D histogram from given xyz files 
    at certain depth D and slab height H.
    NOW WORKS ONLY FOR WATER BEADS.
    """
    Nf = len(dumpfiles)
    res = np.zeros((nbins, nbins))
    axis = set([0, 1, 2]).difference(plane)
    axis = list(axis)[0]    # normal axis to the given plane
    print("Normal axis:", axis)
    for dumpfile in dumpfiles:
        A = read_outfile(dumpfile)
        A = A[(A[:, axis+1] > D-H/2) & (A[:, axis+1] < D+H/2)]
        beads3 = A[A[:, 0] == 3][:, 1:]  # bead C, 3 H2O molecules
        beads4 = A[A[:, 0] == 4][:, 1:]  # bead W, 6 H2O molecules
        profile3, _, _ = np.histogram2d(beads3[:, plane[0]], \
                beads3[:, plane[1]], bins=nbins)
        profile4, _, _ = np.histogram2d(beads4[:, plane[0]], \
                beads4[:, plane[1]], bins=nbins)
        res += 0.5 * profile3 / float(Nf)
        res += 1.0 * profile4 / float(Nf)
    return res


if __name__ == "__main__":
    args = docopt(__doc__)
    L = float(args["--L"])
    axes = {"x": 0, "y": 1, "z": 2}
    planes = {"xy": (0, 1), "yz": (1, 2), "xz": (0, 2)}
    subst_map = {"water": "W", "sulfonic": "S", "backbone": "B", "el": "E"}
    dumpfiles = glob.glob(args["<files>"])
    Nf = len(dumpfiles)
    nbins = int(args["--nbins"])
    subst = args["--subst"]
    if subst not in subst_map.keys():
        sys.exit("Allowed substances: 'water', 'sulfonic', 'backbone', 'el'.")
    if not dumpfiles:
        sys.exit("No files captured, aborting.")

    print("===== Water profile =====")
    print("Substance: %s | Bins: %i | Box size: %.2f | xyz files: %i" % \
            (subst, nbins, L, Nf))
    

    if args["1d"]:
        axis = args["<ax>"]
        if axis not in ["x", "y", "z"]:
            sys.exit("Choose axis from 'x', 'y', 'z'.")
        axis = axes[axis]
        bins = np.linspace(0, L, nbins+1)
        r = bins[:-1] + np.diff(bins)/2.0
        dr = bins[1] - bins[0]

        ti = time.time()
        profile, bins = create_1d_profile(dumpfiles, axis, subst, bins)
        profile /= (dr * L**2)
        tf = time.time()

        outname = "profile_1d_%s.out" % subst_map[subst]
        np.savetxt(outname, np.c_[r, profile])
        print("Done. Time: %.2f s." % (tf - ti))
        print("Array saved in %s." % outname)


    if args["2d"]:
        plane = args["<pl>"]
        if plane not in ["xy", "yz", "xz"]:
            sys.exit("Choose plane from 'xy', 'yz', 'xz'.")
        plane = planes[plane]
        D = float(args["--depth"]) * L
        if D < 0.0 or D > L:
            sys.exit("Choose depth between 0 and 1.")
        H = float(args["--thick"])
        print("2D profile | Plane: %s | Depth: %.1f | Thickness: %.1f") % \
             (args["<plane>"], D, H)

        profile = create_2d_profile(dumpfiles, plane, nbins, D, H)
        outname = "profile_2d.out"
        np.savetxt(outname, profile)
        print("Matrix saved in", outname)

        plt.imshow(profile, cmap="spectral")
        plt.axis("off")
        plotname = "profile_2d_%s.png" % (args["<plane>"],)
        plt.savefig(plotname)
        print("Plot saved in %s.", plotname)


