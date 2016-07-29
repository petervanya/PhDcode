#!/usr/bin/env python
"""Usage:
    gen_binmixt.py <N> <f> [--AB <AB>] [--L <L>] [--m <m>] [--sigma <c>] [--save <fname>]

[AD HOC] create randomly N atoms of 2 types in a box of size 10x10x10
for LAMMPS simulation, atom_style atomic

Arguments:
    <N>              Polymerisation (# of monomers)
    <f>              Fraction of A monomers

Options:
    --AB <AB>        Interaction between different particle types [default: 1.0]
    --L <L>          Box size [default: 10]
    --m <m>          Mass [default: 1.0]
    --sigma <c>      Cutoff for LJ potential [default: 3.0]
    --save <fname>   Save into file

pv278@cam.ac.uk, 21/10/15
"""
import numpy as np
import sys
from docopt import docopt

def get_pos(N, f, L):
    """Return position matrix and type vector of particles"""
    pos = np.random.rand(N, 3)*L
    types = [1]*int(f*N) + [2]*int((1-f)*N)
    return types, pos


def get_header(N, L):
    """Generate LAMMPS header"""
    s = "#blabla\n"
    s += str(N) + " atoms\n"
    s += "2 atom types\n"
    s += "\n"
    s += "0.0 " + str(L) + " xlo xhi\n"
    s += "0.0 " + str(L) + " ylo yhi\n"
    s += "0.0 " + str(L) + " zlo zhi\n\n"
    return s


def get_masses(m1, m2):
    s = "Masses\n\n"
    s += "1 " + str(m1) + "\n"
    s += "2 " + str(m2) + "\n\n"
    return s


def get_pair_coeffs(AA, BB, AB, sigma):
    s = "PairIJ Coeffs # lj/cut\n\n"
    s += "1 1 " + str(AA) + " " + str(sigma) + "\n"
    s += "2 2 " + str(BB) + " " + str(sigma) + "\n"
    s += "1 2 " + str(AB) + " " + str(sigma) + "\n\n"
    return s


def atoms_to_str(types, pos):
    """Convert atomic matrix to string"""
    M = len(pos)
    s = "Atoms\n\n"
    for i in range(M):
        s += "%i\t%i\t%f\t%f\t%f\n" % (i+1, types[i], pos[i, 0], pos[i, 1], pos[i, 2])
    s += "\n"
    return s


if __name__ == "__main__":
    args = docopt(__doc__)
    N = int(args["<N>"])
    f = float(args["<f>"])
    AB = float(args["--AB"])
    L = float(args["--L"])
    m = float(args["--m"])
    sigma = float(args["--sigma"])

    types, pos = get_pos(N, f, L)
    header = get_header(N, L)
    final_string = header + \
                   get_masses(m, m) + \
                   get_pair_coeffs(1.0, 1.0, AB, sigma) + \
                   atoms_to_str(types, pos)

    if args["--save"]:
        fname = args["--save"]
        open(fname, "w").write(final_string)
        print("Data file written in", fname)
    else:
        print(final_string)

