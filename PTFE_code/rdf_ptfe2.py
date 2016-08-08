#!/usr/bin/env python
"""Usage:
    rdf_ptfe2.py <fnames> (--bead <bt> | water) 
                [--L <L> --bins <nbins>]
    rdf_ptfe2.py test

Read xyz files and compute radial distribution function for any 
DPD bead type using Fortran module.
Special treatment of beads containing water, as parametrised by Wu, EES, 2008:
* C: 3 molecules (SO3+3H2O, bead 3)
* W: 6 molecules (bead 4)
TODO
====
Implement RDF of two types.

Arguments:
    --beadtype <bt>  Number from 1 to num. bead types, typically 4

Options:
    --bins <nbins>   Number of bins [default: 500]
    --L <L>          Box size [default: 40.0]

pv278@cam.ac.uk, 05/08/16
"""
import numpy as np
from math import pi
import glob, sys
from docopt import docopt
import lmp_lib as ll
from Fcore.f_rdf import f_rdf      # Fortran module


class mydict(dict):
    """A container for all the system constants"""
    def __getattr__(self, key):
        return self[key]


def guess_box_size(xyz):
    """Infer box size from xyz matrix"""
    return np.round(np.max(xyz[:, 1] - np.min(xyz[:, 1]), 0))


def rdf_water(outfiles, sp):
    """Rdf for water beads from given xyz frames"""
    rdf = np.zeros(len(sp.bins) - 1)
    L = sp.cell[0, 0]
    r = sp.bins[:-1] + np.diff(sp.bins) / 2.0
    dr = r[1] - r[0]

    A = ll.read_xyzfile(outfiles[0])
    NC = len(A[A[:, 0] == 3])
    NW = len(A[A[:, 0] == 4])
    Norm_C = NC * (NC - 1) / 2
    Norm_W = NW * (NW - 1) / 2
    Norm_CW = NC * NW
    mem = NC * NW * 8 / (1024)**3

    print("Computing water. C beads: %i | W beads: %i | Memory: %.1f GB" %\
          (NC, NW, mem))

    for outfile in outfiles:
        A = ll.read_xyzfile(outfile)
        xyz_C = A[A[:, 0] == 3][:, 1:]
        xyz_W = A[A[:, 0] == 4][:, 1:]

        print("Calculating rdf... ", end="", flush=True)
        dist_vec_C = f_rdf.dist_vec(xyz_C, sp.cell)
        rdf_raw_C, _ = np.histogram(dist_vec_C, sp.bins)
        del dist_vec_C
        print("Done.\nCalculating rdf... ", end="", flush=True)
        dist_vec_W = f_rdf.dist_vec(xyz_W, sp.cell)
        rdf_raw_W, _ = np.histogram(dist_vec_W, sp.bins)
        del dist_vec_W
        print("Done.\nCalculating rdf... ", end="", flush=True)
        dist_vec_CW = f_rdf.dist_vec_2mat(xyz_C, xyz_W, sp.cell)
        rdf_raw_CW, _ = np.histogram(dist_vec_CW, sp.bins)
        del dist_vec_CW
        print(outfile, "done.")

        rdf_raw = rdf_raw_C * 3**2 / Norm_C + \
                  rdf_raw_W * 6**2 / Norm_W + \
                  rdf_raw_CW * 3*6 / Norm_CW
        rdf += rdf_raw / (4*pi*r**2 * dr) * L**3 / (3**2 + 6**2 + 3*6)
    return rdf


def rdf_1type(outfiles, sp):
    """Rdf for given bead type from all the available xyz frames"""
    rdf = np.zeros(len(sp.bins) - 1)
    L = sp.cell[0, 0]
    r = sp.bins[:-1] + np.diff(sp.bins) / 2.0
    dr = r[1] - r[0]

    A = ll.read_xyzfile(outfiles[0])
    N1 = len(A[A[:, 0] == sp.beadtype])
    Norm = N1 * (N1-1) / 2

    for outfile in outfiles:
        A = ll.read_xyzfile(outfile)
        xyz = A[A[:, 0] == sp.beadtype][:, 1:]
        print("Calculating rdf... ", end="", flush=True)
        dist_vec = f_rdf.dist_vec(xyz, cell)
        print("Done.\nBinning dist vec... ", end="", flush=True)
        rdf_raw, _ = np.histogram(dist_vec, sp.bins)
        print(outfile, "done.")
        rdf += rdf_raw / (4*pi*r**2 * dr) * L**3 / Norm
    return rdf / sp.Nf


def test_me():
    np.random.seed(1234)
    N = 1000
    L = 1.0
    cell = L * np.eye(3)
    xyz = L * np.random.rand(N, 3)

    Nb = 100
    bins = np.linspace(0, L/2, Nb)
    r = bins[:-1] + np.diff(bins) / 2.0
    dr = r[1] - r[0]
    
    print("Testing rdf. L: %.2f | N: %i | Nbins: %i" % (L, N, Nb))
    dist_vec = f_rdf.dist_vec(xyz, cell)
    rdf, _ = np.histogram(dist_vec, bins)
    rdf = rdf / (4*pi*r**2 * dr) * L**3 / (N*(N-1)/2)
    np.savetxt("test_rdf.out", np.vstack((r, rdf)).T)
    print("Testing rdf of random box saved in test_rdf.out")


if __name__ == "__main__":
    args = docopt(__doc__)
    if args["test"]:
        print("Testing rdf")
        test_me()
        sys.exit()
    
    outfiles = glob.glob(args["<fnames>"])
    if len(outfiles) == 0:
        raise ValueError("No xyz files captured.")
    Nf = len(outfiles)
    N = int(open(outfiles[0], "r").readline())
    Nbins = float(args["--bins"])
    A = ll.read_xyzfile(outfiles[0])
    Nbt = len(set(A[:, 0]))

    if args["--L"]:
        L = float(args["--L"])
    else:
        xyz = ll.read_xyzfile(dumpfiles[0])
        L = guess_box_size(xyz)

    cell = L * np.eye(3)
    bins = np.linspace(0, L/2, Nbins+1)
    r = bins[:-1] + np.diff(bins)/2.0
    
    print("===== Calculating PTFE rdf ====")
    print("Box size: %.2f | Beads: %i | Bins: %i | Frames: %i" %\
         (L, N, Nbins, Nf))


    if args["water"]:
        sp = mydict(N=N, Nf=Nf, cell=cell, bins=bins)
        rdf = rdf_water(outfiles, sp)
        fname = "rdf_water.out"
    elif args["<bt>"]:
        beadtype = int(args["<bt>"])
        if beadtype not in range(1, Nbt+1):
            sys.exit("Unknown bead type, choose 'water' or from 1...%i" % Nbt)
        sp = mydict(N=N, Nf=Nf, cell=cell, bins=bins, beadtype=beadtype)
        rdf = rdf_1type(outfiles, sp)
        fname = "rdf_%i.out" % beadtype
    else:
        sys.exit("Unknown bead type, choose 'water' or from 1...%i" % Nbt)

    np.savetxt(fname, np.vstack((r, rdf)).T)
    print("rdf saved in %s." % fname)


