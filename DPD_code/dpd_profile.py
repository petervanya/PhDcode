#!/usr/bin/env python
"""Usage:
    dpd_profile.py 1d <ax> <frames> --bt <bt> [--L <L> --nbins <nb>]
    dpd_profile.py 2d <ax> <frames> --bt <bt> [--d <d> --h <h> --L <L> --nbins <nb>]

Create 1d or 2d profile of DPD beads.

Arguments:
    <frames>       xyz frames
    <ax>           Axis to profile along ("x", "y", "z")
    <pl>           Plane to profile on ("xy", "yz", "xz")
    --bt <bt>      Bead type

Options:
    --d <d>        Relative depth at which to profile the plane [default: 0.5]
    --h <h>        Relative thickness of the profile [default: 0.1]
    --L <L>        Box size [default: 40]
    --nbins <nb>   Number of bins for the histogram [default: 100]

pv278@cam.ac.uk, 22/02/17
"""
import numpy as np
import matplotlib.pyplot as plt
import glob, sys, time
from dlms_lib import read_xyzfile
from docopt import docopt


def guess_box_size(xyz):
    """Infer box size from xyz matrix"""
    return np.array([np.round(max(xyz[:, i]) - min(xyz[:, i]), 0) \
            for i in range(1, 4)])


def create_1d_profile(frames, bt, ax, bins):
    """Function with pre-set bins and fixed water plot
    * frames: list of frame names
    * ax: x, y, or z
    * bt: bead type
    * bins: vector of bins
    """
    Nf = len(frames)
    res = np.zeros(len(bins)-1)
    for frame in frames:
        A = read_xyzfile(frame)
        beads = A[A[:, 0] == bt][:, 1:]
        profile, _ = np.histogram(beads[:, ax], bins=bins)
        res += profile / Nf
    return res


def create_2d_profile(frames, bt, ax, nbins, D, H):
    """2D histogram from given xyz files 
    at certain depth D and slab height H.
    """
    Nf = len(frames)
    res = np.zeros((nbins, nbins))
    plane = set([0, 1, 2]).difference(ax)
    for frame in frames:
        A = read_xyzfile(frame)
        A = A[(A[:, ax+1] > D-H/2) & (A[:, ax+1] < D+H/2)]
        beads = A[A[:, 0] == bt][:, 1:]
        profile, _, _ = np.histogram2d(beads[:, plane[0]], \
                beads[:, plane[1]], bins=nbins)
        res += profile / float(Nf)
    return res


if __name__ == "__main__":
    args = docopt(__doc__)
    frames = glob.glob(args["<frames>"])
    Nf = len(frames)
    if Nf == 0:
        sys.exit("No frames captured.")
    nbins = int(args["--nbins"])
    bt = int(args["--bt"])
    xyz = read_xyzfile(frames[0])
    xyz_types = set(xyz[:, 0].astype(int))
    if bt not in xyz_types:
        sys.exit("Requested bead type not present in the frames.")
#    L = float(args["--L"])
    L = guess_box_size(xyz)

    ax_raw = args["<ax>"]
    axes = {"x": 0, "y": 1, "z": 2}
    if ax_raw not in axes.keys():
        sys.exit("Choose axis from 'x', 'y', 'z'.")
    ax = axes[ax_raw]

    print("===== DPD profile =====")
    print("Bead type: %i Box: %s" % (bt, str(L)))
    print("Bins: %i | xyz frames: %i" % (nbins, Nf))
    
    if args["1d"]:
        bins = np.linspace(0, L[ax], nbins+1)
        dr = bins[1] - bins[0]
        r = dr / 2 + bins[:-1]

        ti = time.time()
        profile = create_1d_profile(frames, bt, ax, bins)
        Lsurf = L[list(set(range(3)).difference([ax]))]
        print(Lsurf, dr, dr / np.prod(Lsurf))
        profile /= (dr * np.prod(Lsurf))
        tf = time.time()

        outname = "profile_1d_%i_%s.out" % (bt, ax_raw)
        np.savetxt(outname, np.vstack((r, profile)).T)
        print("Time: %.2f s." % (tf - ti))
        print("Profile saved in %s." % outname)


    if args["2d"]:
        D = float(args["--d"])
        H = float(args["--t"])
        if D < 0.0 or D > 1.0 or H < 0.0 or H > 1.0:
            sys.exit("Choose depth between 0 and 1.")
        print("2D profile | Plane: %s | Depth: %.1f | Thickness: %.1f") % \
             (args["<plane>"], D, H)
        D *= L[ax]
        H *= L[ax]

        profile = create_2d_profile(frames, bt, ax, nbins, D, H)
        outname = "profile_2d.out"
        np.savetxt(outname, profile)
        print("Profile saved in %s." % outname)

        plt.imshow(profile, cmap="spectral")
        plt.axis("off")
        plotname = "profile_2d_%i_%s.png" % (bt, ax_raw)
        plt.savefig(plotname)
        print("Plot saved in %s.", plotname)


