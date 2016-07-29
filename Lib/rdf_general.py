#!/usr/bin/env python
"""Usage:
    rdf_general.py <files> [--cutoff <rc> --L <L> --types <n> --bins <nbins>] 

Read xyz files and compute radial distribution function for a given atom type.
Correct implementation of PBCs.
Two options:
* Compute all the pairs
* Consider only pairs up to cutoff <rc> to save memory
  TODO: Cross-distribution does not normalise correctly

Arguments:
    <files>             Regex match xyz files

Options:
    --types <n>         Atom name in number [default: 1]
    --cutoff <rc>       Consider pairs up to <rc> (units of xyz files)
    --bins <nbins>      Number of bins [default: 500]
    --L <L>             Box size (units of xyz files)

pv278@cam.ac.uk, 28/06/16
"""
import numpy as np
from numpy import pi
import glob, sys
from docopt import docopt
from Fcore.f_rdf import f_rdf      # Fortran module
import lmp_lib as ll


class mydict(dict):
    """A container for all the system constants"""
    def __getattr__(self, key):
        return self[key]


def rdf_1type(dumpfiles, sp):
    """Construct an rdf for water beads from given xyz frames.
    Compute radial dist'n fcn from the xyz frame using 
    Fortran routine for pair distances and numpy binning
    for one particle type.
    Input:
    * dumpflies: list of xyz frame names
    * sp: system params, containing:
        * atom_type: integer
        * rc: cutoff for atom neighbourhood
        * cell: (3, 3) matrix
        * bins: vector
        * bc: Bulgarian const to guess the distance vector size
    """
    rdf = np.zeros(len(sp.bins)-1)
    L = sp.cell[0, 0]
    r = sp.bins[:-1] + np.diff(sp.bins)/2.0
    dr = r[1] - r[0]

    A = ll.read_xyzfile(dumpfiles[0])
    N1 = len(A[A[:, 0] == sp.atom_types[0]])
    norm = N1*(N1-1)/2


    for dump in dumpfiles:
        A = ll.read_xyzfile(dump)
        xyz = A[A[:, 0] == sp.atom_types[0]][:, 1:]
        print("Calculating rdf...")

        if sp.use_cutoff:
            nn = int(2* N1**2 * (sp.rc/L)**3 * sp.bc)  # guess of dist vec length
            print("Guessed distance vector size: %i | BG: %.2f" % (nn, sp.bc))
            dist_vec = f_rdf.dist_vec_cut(xyz, sp.rc, L, sp.cell, nn)
            dist_vec = dist_vec[dist_vec != 0.0]
        else:
            dist_vec = f_rdf.dist_vec(xyz, cell)

        rdf_raw, _ = np.histogram(dist_vec, sp.bins)
        rdf += rdf_raw / (4*pi*r**2 * dr) * L**3 / norm
        print("Done: %s" % dump)
    return rdf / sp.Nd


def rdf_2types(dumpfiles, sp):
    """Construct an rdf for water beads from given xyz frames
    Compute radial dist'n fcn from the xyz frame using 
    Fortran routine for pair distances and numpy binning
    for one particle type.
    Similar to rdf_1type."""
    rdf = np.zeros(len(bins)-1)
    L = cell[0, 0]
    r = sp.bins[:-1] + np.diff(sp.bins)/2.0
    dr = r[1] - r[0]
    A = ll.read_xyzfile(dumpfiles[0])
    N1 = len(A[A[:, 0] == sp.atom_types[0]][:, 1:])
    N2 = len(A[A[:, 0] == sp.atom_types[1]][:, 1:])
    norm = N1*N2

    for dump in dumpfiles:
        A = ll.read_xyzfile(dump)
        xyz1 = A[A[:, 0] == sp.atom_types[0]][:, 1:]
        xyz2 = A[A[:, 0] == sp.atom_types[1]][:, 1:]
        print("Atoms: %i %i | Calculating rdf..." % (N1, N2))

        if sp.use_cutoff:
            nn = int(2*max(N1, N2)**2 * (sp.rc/L)**3 * sp.bc)
            print("Guessed distance vector size: %i | BG: %.2f" % (nn, sp.bc))
            dist_vec = f_rdf.dist_vec_cut_2mat(xyz1, xyz2, sp.rc, L, sp.cell, nn)
            dist_vec = dist_vec[dist_vec != 0.0]
        else:
            dist_vec = f_rdf.dist_vec_2mat(xyz1, xyz2, cell)

        rdf_raw, _ = np.histogram(dist_vec, sp.bins)
        rdf += rdf_raw /(4*pi*r**2 * dr) * L**3 / norm
        print("Done: %s" % dump)
    return rdf / sp.Nd


def guess_box_size(xyz):
    """Infer box size from xyz matrix"""
    return np.round(np.max(xyz[:, 1] - np.min(xyz[:, 1]), 0))


if __name__ == "__main__":
    args = docopt(__doc__)
#    print(args)
    dumpfiles = glob.glob(args["<files>"])
    if len(dumpfiles) == 0:
        print("No xyz files captured.")
        sys.exit()
    Nd = len(dumpfiles)
    N = int(open(dumpfiles[0], "r").readline())
    Nbins = int(args["--bins"])

    if args["--L"]:
        L = float(args["--L"])
    else:
        xyz = ll.read_xyzfile(dumpfiles[0])
        L = guess_box_size(xyz)

    use_cutoff = False
    rc = -1.0
    if args["--cutoff"]:
        rc = float(args["--cutoff"])
        use_cutoff = True
    bins = np.linspace(0, L/2, Nbins+1)
    r = bins[:-1] + np.diff(bins)/2.0
    
    atom_types = [int(i) for i in args["--types"].split()]
    if len(atom_types) > 2:
        print("Only two atom types allowed for rdf calculation.")
        sys.exit()
    xyz_types = set(ll.read_xyzfile(dumpfiles[0])[:, 0])
    if not set(atom_types).issubset(xyz_types):
        print("Requested atom types not present in the xyz files.")
        sys.exit()

    cell = L*np.eye(3)
    sp = mydict(N=N, Nd=Nd, cell=cell, bins=bins, atom_types=atom_types, \
                bc=1.3, rc=rc, use_cutoff=use_cutoff)
    
    print("===== Calculating rdf =====")
    print("Atoms: %s | Bins: %i | Box: %.1f | rc %.2f | xyz files: %i" % \
         (repr(atom_types), Nbins, L, rc, Nd))

    if len(atom_types) == 1:
        vals = rdf_1type(dumpfiles, sp)
        fname = "rdf_%i.out" % atom_types[0]
    elif len(atom_types) == 2:
        vals = rdf_2types(dumpfiles, sp)
        fname = "rdf_%i_%i.out" % tuple(atom_types)

    np.savetxt(fname, np.vstack((r, vals)).T, fmt="%.4f")
    print("rdf saved in", fname)


