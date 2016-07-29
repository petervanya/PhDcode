#!/usr/bin/env python
"""Usage:
    dlms_diblock_input.py [--N <N> --f <f> --L <L> --rho <rho>]

Generate diblock copolymer (from A and B particles)
initial configuration and bonds file.

Options:
    --N <N>        Polymerisation [default: 10]
    --f <f>        Fraction of A partcicles [default: 0.5]
    --L <L>        Box size [default: 10.0]
    --rho <rho>    DPD density [default: 3.0]

02/06/16
"""
import numpy as np
from numpy import pi, cos, sin
from docopt import docopt
import sys


def grow_polymer(L, f, n, Nc, mu=1.0, sigma=0.1):
    """Generate coordinates of matrix polymer chains (taken from gen_pmma.py)
    return (n*Nc, 3) xyz matrix
    Input:
    * L: cubic box size
    * n: polymerisation
    * Nc: number of chains
    * mu: mean of bead distance
    * sigma: deviation of bead distance"""
    xyz = np.zeros((Nc*n, 3))
    for i in range(Nc):
        xyz[i*n : (i+1)*n] = grow_one_chain(L, n, Nc, mu, sigma)
    return xyz


def grow_one_chain(L, n, Nc, mu, sigma):
    """Return (n, 3) xyz matrix of one chain"""
    xyz = np.zeros((n, 3))
    xyz[0] = np.random.rand(3)*L
    for i in range(1, n):
        theta = np.random.rand()*pi
        phi = np.random.rand()*2*pi
        r = mu #+ np.random.randn()*L*sigma  # STOPPED USING randn
        new_bead_pos = [r*cos(theta), r*sin(theta)*cos(phi), r*sin(theta)*sin(phi)]
        xyz[i] = xyz[i-1] + new_bead_pos
    return xyz


def gen_bonds(n, Nc):
    """Return (n*Nc, 2) matrix, columns: [atom1, atom2]
    Input:
    * n: polymerisation
    * Nc: number of chains"""
    bond_mat = np.zeros(((n-1)*Nc, 2), dtype=int)
    for i in range(Nc):
        one_chain_bonds = np.hstack(( np.matrix( np.arange(n*i+1, n*(i+1)) ).T,\
                                      np.matrix( np.arange(n*i+2, n*(i+1)+1) ).T ))
        bond_mat[i*(n-1) : (i+1)*(n-1)] = one_chain_bonds
    return bond_mat


def save_config(fname, names, xyz):
    """Save positions into file"""
    N = len(xyz)

    imcon = 0   # include box coordinates
    conf_str = "pokus\n" + "0\t%i\n" % imcon
    if imcon == 1:
        box_mat = L*np.eye(3)
        for i in range(len(box_mat)):
            conf_str += "%f\t%f\t%f\n" % (box_mat[i, 0], box_mat[i, 1], box_mat[i, 2])

    for i in range(N):
        conf_str += "%s        %i\n" % (names[i], i+1)
        conf_str += "    %.10f    %.10f    %.10f\n" % (xyz[i, 0], xyz[i, 1], xyz[i, 2])

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
    N = int(args["--N"])
    f = float(args["--f"])
    L = float(args["--L"])
    rho = float(args["--rho"])
    if f < 0.0 or f > 1.0:
        print("f should be between 0 and 1.")
        sys.exit()
    print("Polymerisation: %i | Fraction: %f" % (N, f))

    Nb = int(rho * L**3)
    Nc = int(Nb/N)
    NA = int(f*N)
    NB = int((1-f)*N)
    k0, r0 = 4.0, 0.1   # spring constants
    
    # ===== positions
    xyz = grow_polymer(L, f, N, Nc)
    names = (["A"]*NA + ["B"]*NB)*Nc
    fname = "CONFIG"
    save_config(fname, names, xyz)

    # ===== bonds
    bond_mat = gen_bonds(N, Nc)
    fname = "bonds.in"
    save_bonds(fname, bond_mat, k0=k0, r0=r0, bond_type="harm")
    

