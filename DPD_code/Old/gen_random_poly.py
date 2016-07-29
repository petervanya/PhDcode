#!/usr/bin/env python
"""Usage:
    gen_random_poly.py (--grow <mu> <sigma> | --mess)
                       <n> <N> <f> [--FAA <AA>] [--FAB <AB>] [--boxsize <L>] 
                       [--seed <s>] [--save <fname>] 

Generate a binary mixture from A and B monomers: AA...AB...B.
Improved version of gen_polymer.py with random polymers.
1. Generate entirely random configurations (mess)
2. Grow the polymer particle by particle (grow)
Uses LAMMPS atom_type molecular

Arguments:
    <n>              Number of molecules
    <N>              Polymerisation (# of monomers)
    <f>              fraction of A monomers
    <mu>             Mean value
    <sigma>          Standard deviation

Options:
    --FAA <AA>       Force between AA (and BB) beads [default: 1.0]
    --FAB <AB>       Force between AB beads [default: 1.0]
    --dist <a>       Initial distance between successive beads
    --save <fname>   Save into file
    --seed <s>       Random seed [default: 123]
    --boxsize <L>    Size of the simulation box L**3 [default: 10]

pv278@cam.ac.uk, 30/10/15
"""
import numpy as np
from docopt import docopt
from xyzlib import Atoms

# ===== physics
def get_pos_mess(N, f=0.5, L=10, mol_num=1):
    """Return position matrix and type vector of particles"""
    #pos = np.random.randn(N, 3)*sigma + mu
    pos = np.random.rand(N, 3)*L
    mol_ids = np.matrix([mol_num]*N).T
    atom_type = np.matrix([1]*int(f*N) + [2]*int((1-f)*N)).T
    return np.hstack((mol_ids, atom_type, pos))


def get_pos_mess2(n, N, f=0.5, L=10):
    """Return position matrix and type vector of particles"""
    pos = np.random.rand(n*N, 3)*L
    atom_type = np.matrix( ([1]*int(f*N) + [2]*int((1-f)*N))*n ).T
    mol_ids = [[i]*N for i in range(1, n+1)]
    mol_ids = np.matrix([item for sublist in mol_ids for item in sublist]).T
    return np.hstack((mol_ids, atom_type, pos))


def get_pos_grow(n, N, f=0.5, L=10, mu=1.0, sigma=0.1):
    """
    For each polymer string
    1. generate first bead
    2. grow the next beads one by one homogenously
    """
    atom_type = np.matrix( ([1]*int(f*N) + [2]*int((1-f)*N))*n ).T
    mol_ids = [[i]*N for i in range(1, n+1)]
    pos = np.zeros((n*N, 3))
    for j in range(0, n):
        pos[j*N] = np.random.randn(3)  # first bead in a string
        for i in range(j*N + 1, (j+1)*N):
            pos[i] = pos[i-1] + np.random.randn(3)*sigma + mu
    mol_ids = [[i]*N for i in range(1, n+1)]
    mol_ids = np.matrix([item for sublist in mol_ids for item in sublist]).T
    return np.hstack((mol_ids, atom_type, pos))


def get_bond_mat(N1, N2):
    """Generate the bond matrix for one molecule:
    [num, type, part1, part2]
    bonds between atom numbers N1, N2 inclusive
    Stiffness constants (type) are the same"""
    N = N2 - N1 + 1  # number of atoms involved
    bond_mat = np.vstack((np.ones(N-1, dtype=int), \
                          np.arange(N1, N2), np.arange(N1+1, N2+1) )).T
    return bond_mat


# ===== printing to string
def get_header(N, Nbonds, L):
    """Generate LAMMPS header"""
    s = "#blabla\n"
    s += str(N) + " atoms\n"
    s += str(Nbonds) + " bonds\n"
    s += "2 atom types\n"
    s += "1 bond types\n"
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


def get_pair_coeffs_LJ(AA, BB, AB, r0):
    s = "PairIJ Coeffs\n\n"
    s += "1 1 " + str(AA) + " " + str(r0) + "\n"
    s += "2 2 " + str(BB) + " " + str(r0) + "\n"
    s += "1 2 " + str(AB) + " " + str(r0) + "\n\n"
    return s


def get_pair_coeffs_DPD(params):   # part1 part2 force gamma cutoff
    s = "PairIJ Coeffs\n\n"
    s += "1 1 " + str(params["1 1"][0]) + " " + str(params["1 1"][1]) + " " +\
         str(params["1 1"][2]) + "\n"
    s += "2 2 " + str(params["2 2"][0]) + " " + str(params["2 2"][1]) + " " +\
         str(params["2 2"][2]) + "\n"
    s += "1 2 " + str(params["1 2"][0]) + " " + str(params["1 2"][1]) + " " +\
         str(params["1 2"][2]) + "\n\n"
#    s += "1 1 " + str(AA) + " 1.0 \n"
#    s += "2 2 " + str(BB) + " 1.0 \n"
#    s += "1 2 " + str(AB) + " 1.0 \n\n"
    return s


def get_bond_coeffs(id_bondcoeffs):
    """Save bond coefficients into a string"""
    s = "Bond Coeffs\n\n"
    for k, v in id_bondcoeffs.iteritems():
        str_v = ""
        for i in len(v):
            str_v += str(i) + "\t"
        s += str(k) + "\t" + str_v + "\n"
    return s


def atoms_to_str(xyz_mat, mol_id=1):
    """Convert atomic matrix to string, atom_type molecular
    atom_mat[:, 0] are atom ids"""
    M = len(xyz_mat)
    s = ""
    for i in range(M):
        s += "%i\t%i\t%i\t%f\t%f\t%f\n" % \
             (i+1, mol_id, xyz_mat[i, 0], xyz_mat[i, 1], xyz_mat[i, 2], xyz_mat[i, 3])
    return s


def atoms_to_str2(atom_mat):
    """Convert atomic matrix to string, atom_type molecular
    xyz_mat[:, 0] are atom ids"""
    M = len(atom_mat)
    s = ""
    for i in range(M):
        s += "%i\t%i\t%i\t%f\t%f\t%f\n" % \
             (i+1, atom_mat[i, 0], atom_mat[i, 1], atom_mat[i, 2], atom_mat[i, 3], atom_mat[i, 4])
    return s


def bonds_to_str(bond_mat):
    """Convert bond_matrix to string"""
    M, N = bond_mat.shape
    s = ""
    for i in range(M):
        s += str(i+1) + "\t"
        for j in range(N):
            s += str(bond_mat[i, j]) + "\t"
        s += "\n"
    return s


if __name__ == "__main__":
    args = docopt(__doc__)
#    print args
    n = int(args["<n>"])
    N = int(args["<N>"])
    f = float(args["<f>"])
    L = float(args["--boxsize"])
    FAA = float(args["--FAA"])
    FAB = float(args["--FAB"])
    fname = args["--save"]
    np.random.seed(int(args["--seed"]))
# DPD PARAMS
    m = 1.0          # PERHAPS SHOULD BE GIVEN AS CMD LN ARG
    gamma = 1.0
    cutoff = 1.0
    pair_ij_coeffs = {"1 1": [FAA, gamma, cutoff],
                      "2 2": [FAA, gamma, cutoff],
                      "1 2": [FAB, gamma, cutoff]}
    
    final_atoms_str = ""  # for all molecules
    final_bonds_str = ""
    if args["--mess"]:
#        atoms_mat = get_pos_mess(N, f, L)
#        final_atoms_str += atoms_to_str(atoms_mat, mol_id)
        atoms_mat = get_pos_mess2(n, N, f, L)
        final_atoms_str += atoms_to_str2(atoms_mat)
    elif args["--grow"]:
        mu = float(args["<mu>"])
        sigma = float(args["<sigma>"])
        atoms_mat = get_pos_grow(n, N, f, L, mu, sigma)
        final_atoms_str += atoms_to_str2(atoms_mat)
    final_atoms_str += "\n"

    for i in range(1, n+1):  # loop through mol_ids
        bonds_mat = get_bond_mat(N*(i-1) + 1, N*i)
        final_bonds_str += bonds_to_str(bonds_mat)

    final_string = get_header(N*n, n*(N-1), L) + \
                   get_masses(m, m) + \
                   get_pair_coeffs_DPD(pair_ij_coeffs) + \
                   "Atoms\n\n" + final_atoms_str + \
                   "Bonds\n\n" + final_bonds_str

    if args["--save"]:
        fname = args["--save"]
        open(fname, "w").write(final_string)
        print "Data file written in", fname
    else:
        print final_string


