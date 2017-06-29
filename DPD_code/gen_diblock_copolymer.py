#!/usr/bin/env python
"""Usage:
    gen_diblock_copolymer.py [--N <N> --f <f> --L <L> --rho <rho> --xyz <xyz>]

Generate a binary mixture from A and B monomers: AA...ABB...B

Options:
    --N <N>        Polymerisation [default: 10]
    --f <f>        Fraction of A partcicles [default: 0.5]
    --L <L>        Box size, 'L' or 'Lx Ly Lz' [default: 10.0]
    --rho <rho>    DPD density [default: 3.0]

pv278@cam.ac.uk, 29/06/16
"""
import numpy as np
from math import sin, cos, acos, pi
from docopt import docopt
import lmp_lib as ll


def grow_polymer(L, f, N, Nc, mu=1.0):
    """Generate coordinates of matrix polymer chains (taken from gen_pmma.py)
    return (n*Nc, 5) matrix, columns: [mol_ids, bead_type, xyz]
    Input:
    * L: cubic box size
    * N: polymerisation
    * Nc: number of chains
    * mu: mean of bead distance"""
    xyz = np.zeros((Nc*N, 5))
    NA, NB = int(N * f), N - int(N * f)
    atom_ids = np.array([1] * NA + [2] * NB)
    for i in range(Nc):
        xyz[i*N : (i+1)*N] = np.c_[(i+1)*np.matrix(np.ones(N)).T, atom_ids, \
                grow_one_chain(L, N, Nc, mu)]
    return xyz


def grow_one_chain(L, n, Nc, mu=1.0):
    """Return (n, 3) xyz matrix of one chain (taken from gen_pmma.py)"""
    xyz = np.zeros((n, 3))
    xyz[0] = np.random.rand(3) * L
    for i in range(1, n):
        th = acos(1 - 2 * np.random.rand())
        phi = 2 * pi * np.random.rand()
        new_pos = [mu * sin(th) * cos(phi), mu * sin(th) * sin(phi), \
                mu * cos(th)]
        xyz[i] = xyz[i-1] + new_pos
    return xyz


def gen_bonds(n, Nc):
    """Create bond matrix (taken from gen_pmma.py)
    return (n*Nc, 3) matrix, columns: [bond_type, atom1, atom2]
    Input:
    * n: polymerisation
    * Nc: number of chains"""
    mat = np.zeros(((n-1)*Nc, 3), dtype=int)
    for i in range(Nc):
        one_chain_bonds = np.c_[np.ones(n-1), \
                np.arange(n*i+1, n*(i+1)), np.arange(n*i+2, n*(i+1)+1)]
        mat[i*(n-1) : (i+1)*(n-1)] = one_chain_bonds
    return mat


if __name__ == "__main__":
    args = docopt(__doc__)
    N = int(args["--N"])
    f = float(args["--f"])
    rho = float(args["--rho"])
    np.random.seed(1234)
    s = args["--L"].split()
    if len(s) == 1:
        L = float(s[0]) * np.ones(3)
    elif len(s) == 3:
        L = np.array(s).astype(float)
    else:
        sys.exit("<L> should have size 1 or 3.")
    if f < 0.0 or f > 1.0:
        sys.exit("f should be between 0 and 1.")

    Nb = int(rho * np.prod(L))
    Nc = int(Nb / N)
    rc = 1.0
    
    print("=== Creating LAMMPS data file for diblock copolymer melt ===")
    print("Set interaction params in the input file")
    print("Box: %s | Rho: %.1f | Chain length: %i | A beads/chain: %i" % \
          (L, rho, N, int(N*f) ))
    
    poly_xyz = grow_polymer(L, f, N, Nc, mu=1.0)
    xyz_str = ll.atoms2str(poly_xyz)
    print(len(poly_xyz), "beads created, density:", len(poly_xyz) / L**3)

    bonds = gen_bonds(N, Nc)
    bonds_str = ll.bonds2str(bonds)
    print(len(bonds), "bonds created")

    final_string = ll.header2str(len(poly_xyz), len(bonds), 2, 1, L[0]) + \
                   ll.mass2str({1: 1.0, 2: 1.0}) + \
                   "\nAtoms\n\n" + xyz_str +  \
                   "\nBonds\n\n" + bonds_str
    
    fname = "diblock.data"
    open(fname, "w").write(final_string)
    print("Data file written in", fname)

    if args["--xyz"]:
        fname = args["--xyz"]
        ll.save_xyzfile(fname, poly_xyz[:, 1:])
        print("xyz file saved in", fname)


