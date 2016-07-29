#!/usr/bin/env python
"""
A collection of functions to manipulate LAMMPS files

pv278@cam.ac.uk, 11/01/16
"""
import numpy as np

# ===== print input
def header2str(N, Nbonds, atomtypes, bondtypes, L):
    """Generate LAMMPS header"""
    s = "#blabla\n"
    s += str(N) + " atoms\n"
    s += str(Nbonds) + " bonds\n"
    s += str(atomtypes) + " atom types\n"
    s += str(bondtypes) + " bond types\n"
    s += "\n"
    s += "0.0 %.1f xlo xhi\n" % L
    s += "0.0 %.1f ylo yhi\n" % L
    s += "0.0 %.1f zlo zhi\n\n" % L
    return s


def mass2str(masses):
    """Print mass dictionary into string for LAMMPS data file"""
    s = "Masses\n\n"
    for k, v in masses.items():
        s += "%i %.2f\n" % (k, v)
    return s + "\n"


def pair_dpd_coeffs2str(coeffs):
    """
    Structure:
    * key: "part1 part2"
    * value: [force, gamma, cutoff]
    """
    s = "PairIJ Coeffs\n\n"
    for k, v in coeffs.items():
        s += "%s %.3f %.2f %.2f\n" % (k, v[0], v[1], v[2])
    return s + "\n"


def bond_coeffs2str(k_ij):
    """Print bond coefficients into string.
    Structure:
    * key: 1..4
    * value [k_ij, r0]
    """
    s = "Bond Coeffs\n\n"
    for k, v in k_ij.items():
        s += "%s %.2f %.2f\n" % (k, v[0], v[1])
    return s + "\n"


def atoms2str(mat, atom_style="molecular"):
    """Convert atomic matrix to str, atom_type molecular
    xyz_mat[:, 0] are atom ids"""
    M = len(mat)
    s = ""
    if atom_style == "molecular":
        for i in range(M):
            s += "%i\t%i\t%i\t%e\t%e\t%e\n" % \
                 (i+1, mat[i, 0], mat[i, 1], mat[i, 2], mat[i, 3], mat[i, 4])
    elif atom_style == "atomic":
        for i in range(M):
            s += "%i\t%i\t%e\t%e\t%e\n" % \
                 (i+1, mat[i, 0], mat[i, 1], mat[i, 2], mat[i, 3])
    return s + "\n"


def bonds2str(bond_mat):
    """Convert bond matrix to string"""
    M, N = bond_mat.shape
    s = ""
    for i in range(M):
        s += str(i+1) + "\t"
        for j in range(N):
            s += str(bond_mat[i, j]) + "\t"
        s += "\n"
    return s + "\n"


def bonds2str2(bond_mat):
    """Take (N, 2) matrix and convert it to LAMMPS string.
    stick single bond id (1) to and bond number to the left.
    Columns: [num, bond_id, bond1, bond2]"""
    M = len(bond_mat)
    s = ""
    for i in range(M):
        s += "%i\t%i\t%i\t%i\n" % (i+1, 1, bond_mat[i, 0], bond_mat[i, 1])
    return s + "\n"


# ===== manipulate output
def read_xyzfile(outfile):
    """Read one xyz outfile into a numpy matrix"""
    A = open(outfile, "r").readlines()[2:]
    A = [line.split() for line in A]
    A = np.array(A, order="F").astype(float)
    return A


def save_xyzfile(fname, mat):
    """Take xyz matrix [ids, x, y, z] and save into fname"""
    N = len(mat)
    with open(fname, "w") as f:
        f.write(str(N) + "\nbla\n")
        for i in range(N):
            f.write("%i\t%f\t%f\t%f\n" % (mat[i, 0], mat[i, 1], mat[i, 2], mat[i, 3]))


