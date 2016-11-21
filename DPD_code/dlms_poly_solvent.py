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

pv278@cam.ac.uk, 18/11/16
"""
import numpy as np
from numpy import pi, cos, sin
import sys
from docopt import docopt
from dlms_lib import save_config, inter2str, species2str2, mol2str


def grow_polymer(L, f, r, Nc, mu=0.8):
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
        xyz[i*r : (i+1)*r] = grow_one_chain(L, r, Nc, mu)
    return xyz


def grow_one_chain(L, r, Nc, mu):
    """Return (n, 3) xyz matrix of one chain"""
    xyz = np.zeros((r, 3))
    xyz[0] = np.random.rand(3) * L
    for i in range(1, r):
        th = np.random.rand() * pi
        ph = np.random.rand() * 2 * pi
        r = mu
        xyz[i] = xyz[i-1] + [r*cos(th), r*sin(th)*cos(ph), r*sin(th)*sin(ph)]
    return xyz


def gen_bonds(n, Nc):
    """Return (n*Nc, 2) matrix, columns: [atom1, atom2]
    Input:
    * n: polymerisation
    * Nc: number of chains"""
    bond_mat = np.zeros(((n-1)*Nc, 2), dtype=int)
    for i in range(Nc):
        one_chain = np.hstack(( np.matrix( np.arange(n*i+1, n*(i+1)) ).T,\
                                np.matrix( np.arange(n*i+2, n*(i+1)+1) ).T ))
        bond_mat[i*(n-1) : (i+1)*(n-1)] = one_chain
    return bond_mat


def gen_bonds2(r):
    """Return (r, 2) matrix, columns: [atom1, atom2]
    Input:
    * r: polymerisation
    * Nc: number of chains"""
    return np.vstack((np.arange(1, r), np.arange(2, r+1))).T


def save_config(fname, names, xyz):
    """Save positions into file"""
    N = len(xyz)

    imcon = 0   # include box coordinates
    conf_str = "bla\n" + "0\t%i\n" % imcon
    if imcon == 1:
        box_mat = L*np.eye(3)
        for i in range(len(box_mat)):
            conf_str += "%f\t%f\t%f\n" % \
                    (box_mat[i, 0], box_mat[i, 1], box_mat[i, 2])

    for i in range(N):
        conf_str += "%s        %i\n" % (names[i], i+1)
        conf_str += "    %.10f    %.10f    %.10f\n" % \
                (xyz[i, 0], xyz[i, 1], xyz[i, 2])

    open(fname, "w").write(conf_str)
    print("Initial configuration written in %s" % fname)


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
    L = float(args["--L"])
    rho = float(args["--rho"])
    aii = float(args["--aii"])
    da = float(args["--da"])
    if f < 0.0 or f > 1.0:
        sys.exit("f should be between 0 and 1.")

    N = int(rho * L**3)
    Nc = int(N * f // r)
    Ns = N - Nc
    k0, r0 = 4.0, 0.1      # spring constants
    print("Fraction of P: %.2f | Num. chains: %i | Polymerisation: %i" % \
            (f, Nc, r))
    
    # ===== positions
    xyzP = grow_polymer(L, f, r, Nc)
    xyzS = np.random.rand(Ns, 3) * L
    xyz = np.vstack((xyzP, xyzS))
    names = ["C"] * r * Nc + ["S"] * Ns
    save_config("CONFIG", names, xyz)

    # ===== generate FIELD file with bead species and interactions
    rc, gamma = 1.0, 4.5
    a_ij = {}
    a_ij["C C"] = [aii, rc, gamma]
    a_ij["S S"] = [aii, rc, gamma]
    a_ij["C S"] = [aii + da, rc, gamma]

    bond_mat = gen_bonds2(r)
    bead_list = ["C"] * r
    mol_str = mol2str("poly", Nc, bead_list, bond_mat, \
                                  bond_type="harm", k0=k0, r0=r0)

    field_string = "bla\n\n" +\
                   species2str2({"C": 0, "S": Ns}) +\
                   inter2str(a_ij) +\
                   "MOLECULES 1\n" + \
                   mol_str + "\n" + \
                   "close\n"
    open("FIELD", "w").write(field_string)
    print("FIELD file saved.")


