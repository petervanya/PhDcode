#!/usr/bin/env python
"""Usage:
    dlms_poly_solvent.py [--r <r> --f <f> --L <L> --rho <rho>]
                         [--aii <aii> --da <da>]

Generate polymer (P beads of chain length r) in solvent (S beads).
Initial CONFIG and FIELD file.

Options:
    --r <r>        Polymerisation [default: 10]
    --f <f>        Fraction of A partcicles [default: 0.5]
    --L <L>        Box size [default: 10.0]
    --rho <rho>    DPD density [default: 3.0]
    --aii <aii>    Default bead repulsion [default: 25]
    --da <da>      Excess repulsion [default: 1]

pv278@cam.ac.uk, 18/11/16, modified 15/02/16
"""
import numpy as np
from math import pi, cos, sin, acos
import sys
from docopt import docopt
from dlms_lib import save_config, inter2str, species2str, mol2str


def grow_polymer(L, f, r, Nc, mu=1.0):
    """Generate coordinates of matrix polymer chains.
    Return (r * Nc, 3) xyz matrix
    Input:
    * L: cubic box size
    * r: polymerisation
    * Nc: number of chains
    * mu: mean bead distance
    * sigma: deviation of bead distance"""
    xyz = np.zeros((Nc*r, 3))
    for i in range(Nc):
        xyz[i*r : (i+1)*r] = grow_one_chain(L, r, mu)
    return xyz


def grow_one_chain(L, r, mu):
    """Return (n, 3) xyz matrix of one chain"""
    xyz = np.zeros((r, 3))
    xyz[0] = np.random.rand(3) * L
    for i in range(1, r):
        th = acos(1 - 2 * np.random.rand())
        ph = 2 * pi * np.random.rand()
        r = mu
        xyz[i] = xyz[i-1] + [r*cos(th), r*sin(th)*cos(ph), r*sin(th)*sin(ph)]
    return xyz


def gen_bonds(r):
    """Return (r, 2) matrix, columns: [atom1, atom2]
    * r: polymerisation"""
    return np.c_[np.arange(1, r), np.arange(2, r+1)]


def save_bonds(fname, bond_mat, k0=4.0, r0=0.1, bond_type="harm"):
    """Save bond matrix into a file.
    This file should be copied into an appropriate place in FIELDS file"""
    bonds_str = ""
    N = len(bond_mat)
    for i in range(N):
        bonds_str += "%s  %i  %i  %f  %f\n" %\
                     (bond_type, bond_mat[i, 0], bond_mat[i, 1], k0, r0)
    open(fname, "w").write(bonds_str)
    print("Bonds saved in %s. Include this matrix in FIELDS file." % fname)


if __name__ == "__main__":
    args = docopt(__doc__)
    r = int(args["--r"])   # polymerisation
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

    N = int(rho * np.prod(L)) // 10 * 10     # round to 10
    box = np.diag(L)
    Nc = int(N * f // r)
    Ns = N - Nc * r

    print("===== DPD polymer-solvent mixture =====")
    print("rho: %.1f | Box: %s" % (rho, str(L)))
    print("Fraction of P: %.2f | Num. chains: %i | Polymerisation: %i" % \
            (f, Nc, r))
    print("Solvent particles: %i" % Ns)
    
    rc, gamma = 1.0, 4.5
    a_ij = {}
    a_ij["C C"] = [aii, rc, gamma]
    a_ij["S S"] = [aii, rc, gamma]
    a_ij["C S"] = [aii + da, rc, gamma]

    # ===== positions
    xyzP = grow_polymer(L, f, r, Nc)
    xyzS = np.random.rand(Ns, 3) * L
    xyz = np.vstack((xyzP, xyzS))
    names = ["C"] * r * Nc + ["S"] * Ns
    save_config("CONFIG", names, xyz, box)

    # ===== bonds
    k0, r0 = 4.0, 0.1      # spring constants
    bond_mat = gen_bonds(r)
    bead_list = ["C"] * r
    mol_str = mol2str("poly", Nc, bead_list, bond_mat, \
            bond_type="harm", k0=k0, r0=r0)

    field_string = "bla\n\n" + \
            species2str({"C": 0, "S": Ns}) + \
            inter2str(a_ij) + \
            "MOLECULES 1\n" + \
            mol_str + "\n" + \
            "close\n"
    open("FIELD", "w").write(field_string)
    print("FIELD file saved.")


