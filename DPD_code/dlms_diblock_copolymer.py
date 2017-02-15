#!/usr/bin/env python
"""Usage:
    dlms_diblock_copolymer.py [--N <N> --f <f> --L <L> --rho <rho>]
                              [--aii <aii> --da <da>]

Generate diblock copolymer (from A and B particles)
initial configuration and bonds file.

Options:
    --N <N>        Polymerisation [default: 10]
    --f <f>        Fraction of A partcicles [default: 0.5]
    --L <L>        Box size, 'L' or 'Lx Ly Lz' [default: 10.0]
    --rho <rho>    DPD density [default: 3.0]
    --aii <aii>    Default bead repulsion [default: 25]
    --da <da>      Excess repulsion [default: 1]

pv278@cam.ac.uk, 15/02/16
"""
import numpy as np
from math import pi, cos, sin, acos
import sys
from dlms_lib import save_config, inter2str, species2str, mol2str
from docopt import docopt


def grow_polymer(L, f, N, Nc, mu=1.0, sigma=0.1):
    """Generate coordinates of matrix polymer chains (taken from gen_pmma.py)
    return (N*Nc, 3) xyz matrix
    Input:
    * L: cubic box size
    * N: polymerisation
    * Nc: number of chains"""
    xyz = np.zeros((Nc*N, 3))
    for i in range(Nc):
        xyz[i*N : (i+1)*N] = grow_one_chain(L, N, mu=mu)
    return xyz


def grow_one_chain(L, N, mu=1.0):
    """Return (N, 3) xyz matrix of one chain"""
    xyz = np.zeros((N, 3))
    xyz[0] = np.random.rand(3) * L
    for i in range(1, N):
        th = acos(1 - 2 * np.random.rand())
        phi = 2 * pi * np.random.rand()
        new_pos = [mu * cos(th), mu * sin(th) * cos(phi), \
                mu * sin(th) * sin(phi)]
        xyz[i] = xyz[i-1] + new_pos
    return xyz


def gen_bonds(N):
    """Return (N, 2) matrix, columns: [atom1, atom2]
    * N: polymerisation"""
    return np.vstack((np.arange(1, N), np.arange(2, N+1))).T


if __name__ == "__main__":
    args = docopt(__doc__)
    N = int(args["--N"])
    f = float(args["--f"])
    rho = float(args["--rho"])
    aii = float(args["--aii"])
    da = float(args["--da"])
    s = args["--L"].split()
    if len(s) == 1:
        L = float(s[0]) * np.ones(3)
    elif len(s) == 3:
        L = np.array(s).astype(float)
    else:
        sys.exit("<L> should have size 1 or 3.")
    if f < 0.0 or f > 1.0:
        sys.exit("f should be between 0 and 1.")

    print("===== DPD diblock copolymer =====")
    print("rho: %.1f | Box: %s" % (rho, str(L)))
    print("Polymerisation: %i | Fraction: %.2f" % (N, f))

    Nb = int(rho * np.prod(L)) // 10 * 10     # round to 10
    box = np.diag(L)
    Nc = int(Nb / N)
    NA = int(f * N)
    NB = int((1 - f) * N)
    
    rc, gamma = 1.0, 4.5
    a_ij = {}
    a_ij["A A"] = [aii, rc, gamma]
    a_ij["B B"] = [aii, rc, gamma]
    a_ij["A B"] = [aii + da, rc, gamma]

    # ===== positions
    xyz = grow_polymer(L, f, N, Nc)
    bead_list = ["A"] * NA + ["B"] * NB
    names = bead_list * Nc
    fname = "CONFIG"
    save_config(fname, names, xyz, box)

    # ===== bonds
    k0, r0 = 4.0, 0.1   # spring constants
    bond_mat = gen_bonds(N)
    mol_str = mol2str("diblock", Nc, bead_list, bond_mat, \
            bond_type="harm", k0=k0, r0=r0)

    field_string = "bla\n\n" + \
            species2str({"A": 0, "B": 0}) +\
            inter2str(a_ij) + \
            "MOLECULES 1\n" + \
            mol_str + "\n" + \
            "close\n"

    fname = "FIELD"
    open(fname, "w").write(field_string)
    print("FIELD file saved.")


