#!/usr/bin/env python
"""
[AD HOC] For Material Studio output files.
Substitute bead names A, B11, C, W, G for numbers 1..5
and convert to metres (now AA).

Usage:
    subst_names.py <input>

02/04/15
"""
import numpy as np
import glob
from docopt import docopt

names_dict = {"A": "1", "B11": "2", "C": 3, "W": 4, "G": 5, "P": 6}


def read_xyzfile(outfile):
    """Read one xyz outfile into a numpy matrix"""
    A = open(outfile, "r").readlines()[2:]
    A = [line.split() for line in A]
    A = np.array(A, order="F")
    return A


def save_xyzfile(fname, mat):
    """Take xyz matrix [ids, x, y, z] and save into fname"""
    N = len(mat)
    with open(fname, "w") as f:
        f.write(str(N) + "\nbla\n")
        for i in range(N):
            f.write("%i\t%e\t%e\t%e\n" % \
                   (mat[i, 0], mat[i, 1], mat[i, 2], mat[i, 3]))


args = docopt(__doc__)
infiles = glob.glob(args["<input>"])
print(infiles)
outfiles = ["dump_" + fn[1:] + ".xyz" for fn in infiles]

for i in range(len(infiles)):
    A = read_xyzfile(infiles[i])
    for ch, d in names_dict.items():
        A[A[:, 0] == ch, 0] = d

    A = A.astype(float)
    A[:, 1:] *= 1e-10
    save_xyzfile(outfiles[i], A)
    print("Converted file saved in", outfiles[i])

