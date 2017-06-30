#!/usr/bin/env python
"""Usage:
    rdf_general.py <frames> --bt <bt> [--rc <rc> --L <L> --bins <nb> --verbose]

Read xyz files and compute radial distribution function for a given atom type.
Correct implementation of PBCs.
Two options:
* Compute all the pairs
* Consider only pairs up to cutoff <rc> to save memory
  TODO: Cross-distribution does not normalise correctly

Arguments:
    <frames>       Regex for xyz frames
    --bt <bt>      Atom name in number

Options:
    --rc <rc>      Consider pairs up to <rc> (units of xyz files)
    --bins <nb>    Number of bins [default: 500]
    --L <L>        Box size (units of xyz files)

pv278@cam.ac.uk, 28/06/16
"""
import numpy as np
from numpy import pi
import glob, sys, time
from docopt import docopt
from Fcore.f_rdf import f_rdf      # Fortran module
from dlms_lib import read_xyzfile


class mydict(dict):
    """A container for all the system constants"""
    def __getattr__(self, key):
        return self[key]


def rdf_1type(frames, sp):
    """Construct an rdf for water beads from given xyz frames.
    Compute radial dist'n fcn from the xyz frame using 
    Fortran routine for pair distances and numpy binning
    for one particle type.
    Input:
    * frames: list of xyz frame names
    * sp: system params, containing:
        * atom_type: integer
        * rc: cutoff for atom neighbourhood
        * box: (3, 3) matrix
        * bins: vector
        * bc: Bulgarian const to guess the distance vector size
    """
    rdf = np.zeros(len(sp.bins)-1)
    L = sp.box[0, 0]
    V = np.prod(np.diag(sp.box))
    r = sp.bins[:-1] + np.diff(sp.bins)/2.0
    dr = r[1] - r[0]

    A = read_xyzfile(frames[0])
    N1 = len(A[A[:, 0] == sp.atom_types[0]])
    norm = N1*(N1-1)/2


    for frame in frames:
        A = read_xyzfile(frame)
        xyz = A[A[:, 0] == sp.atom_types[0]][:, 1:]
        if sp.verbose: print("Calculating rdf... ", end="")

        if sp.use_cutoff:
            nn = int(2* N1**2 * (sp.rc/L)**3 * sp.bc)  # guess dist vec size
            if sp.verbose:
                print("Guessed dist. vec. size: %i | BG: %.2f" % (nn, sp.bc))
            dist_vec = f_rdf.dist_vec_cut(xyz, sp.rc, L, sp.box, nn)
            dist_vec = dist_vec[dist_vec != 0.0]
        else:
            dist_vec = f_rdf.dist_vec(xyz, sp.box)

        rdf_raw, _ = np.histogram(dist_vec, sp.bins)
        rdf += rdf_raw / (4 * pi * r**2 * dr) * V / norm
        if sp.verbose: print("Done: %s" % frame)
    return rdf / sp.Nf


def rdf_2types(frames, sp):
    """Construct an rdf for water beads from given xyz frames
    Compute radial dist'n fcn from the xyz frame using 
    Fortran routine for pair distances and numpy binning
    for one particle type.
    Similar to rdf_1type."""
    rdf = np.zeros(len(bins)-1)
    L = sp.box[0, 0]
    V = np.prod(np.diag(sp.box))
    r = sp.bins[:-1] + np.diff(sp.bins)/2.0
    dr = r[1] - r[0]
    A = read_xyzfile(frames[0])
    N1 = len(A[A[:, 0] == sp.atom_types[0]][:, 1:])
    N2 = len(A[A[:, 0] == sp.atom_types[1]][:, 1:])
    norm = N1 * N2

    for frame in frames:
        A = read_xyzfile(frame)
        xyz1 = A[A[:, 0] == sp.atom_types[0]][:, 1:]
        xyz2 = A[A[:, 0] == sp.atom_types[1]][:, 1:]
        if sp.verbose:
            print("Atoms: %i %i | Calculating rdf... " % (N1, N2), end="")

        if sp.use_cutoff:
            nn = int(2*max(N1, N2)**2 * (sp.rc/L)**3 * sp.bc)
            print("Guessed distance vector size: %i | BG: %.2f" % (nn, sp.bc))
            dist_vec = f_rdf.dist_vec_cut_2mat(xyz1, xyz2, sp.rc, L, sp.box, nn)
            dist_vec = dist_vec[dist_vec != 0.0]
        else:
            dist_vec = f_rdf.dist_vec_2mat(xyz1, xyz2, box)

        rdf_raw, _ = np.histogram(dist_vec, sp.bins)
        rdf += rdf_raw /(4 * pi * r**2 * dr) * V / norm
        if sp.verbose: print("Done: %s" % frame)
    return rdf / sp.Nf


def guess_box_size(A):
    """Infer box size from xyz matrix"""
    return np.array([np.round(np.max(A[:, i]) - np.min(A[:, i]), 2) \
            for i in range(1, 4)])


def guess_box_size2(frames):
    """Infer the box size from the xyz frame"""
    Ntf = len(frames) // 50 + 1
    Ltrial = np.zeros((Ntf, 3))
    for i in range(Ntf):
        A = read_xyzfile(frames[i])
        Ltrial[i] = np.array([np.round(np.max(A[:, i]) - np.min(A[:, i]), 2) \
                for i in range(1, 4)])
    return np.max(Ltrial, 0)


if __name__ == "__main__":
    args = docopt(__doc__)
    frames = glob.glob(args["<frames>"])
    if len(frames) == 0:
        sys.exit("No xyz frames captured.")
    Nf = len(frames)
    N = int(open(frames[0], "r").readline())
    Nbins = int(args["--bins"])

    if args["--L"]:
        s = args["--L"].split()
        if len(s) == 1:
            Ls = float(s[0]) * np.ones(3)
        elif len(s) == 3:
            Ls = np.array(s).astype(float)
        else:
            sys.exit("L should have one or three elements.")
    else:
        Ls = guess_box_size2(frames)

    use_cutoff = False
    rc = -1.0
    if args["--rc"]:
        rc = float(args["--rc"])
        use_cutoff = True
    bins = np.linspace(0, np.min(Ls) / 2, Nbins+1)
    dr = bins[1] - bins[0]
    r = dr / 2 + bins[:-1]
    
    atom_types = [int(i) for i in args["--bt"].split()]
    if len(atom_types) > 2:
        sys.exit("Only two atom types allowed for rdf calculation.")
    xyz_types = set(read_xyzfile(frames[0])[:, 0])
    if not set(atom_types).issubset(xyz_types):
        sys.exit("Requested atom types not present in the xyz files.")

    box = np.diag(Ls)
    sp = mydict(N=N, Nf=Nf, box=box, bins=bins, atom_types=atom_types, \
                bc=1.3, rc=rc, use_cutoff=use_cutoff, \
                verbose=args["--verbose"])
    
    print("===== Calculating rdf =====")
    print("Atoms: %s | Bins: %i | rc %.2f | xyz frames: %i" % \
         (repr(atom_types), Nbins, rc, Nf))
    print("Box: %s" % Ls)

    ti = time.time()
    if len(atom_types) == 1:
        vals = rdf_1type(frames, sp)
        fname = "rdf_%i.out" % atom_types[0]
    elif len(atom_types) == 2:
        vals = rdf_2types(frames, sp)
        fname = "rdf_%i_%i.out" % tuple(atom_types)
    tf = time.time()
    print("Time: %.2f s." % (tf - ti))

    np.savetxt(fname, np.c_[r, vals], fmt="%.4f")
    print("RDF saved in %s." % fname)


