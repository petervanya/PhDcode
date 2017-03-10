#!/usr/bin/env python
"""Usage: dpd_pressure.py <frame> [--lc]

Compute the virial pressure.

Options:
    --lc      Use link cells

pv278@cam.ac.uk, 06/03/17
"""
import numpy as np
from numpy.linalg import norm, pinv
from math import sqrt
from numba import jit, float64, int64
import sys, time
from docopt import docopt
from dlms_lib import read_xyzfile
from link_cells import LinkCells


def parse_field():
    """Extract coefficients from the field file.
    Assume rc = 1.0"""
    try:
        f = open("FIELD", "r").readlines()
    except FileNotFoundError:
        sys.exit("FIELD file not found in the current dir.")
    for line in range(len(f)):
        if "SPECIES" in f[line].upper():
            Nb = int(f[line].split()[-1])
            sp_lines = f[line+1 : line+Nb+1]
        if "INTERACTIONS" in f[line].upper():
            Ni = int(f[line].split()[-1])
            int_lines = f[line+1 : line+Ni+1]

    sp = {}       # species
    sp_num = {}   # species dict with numbers
    for i in range(Nb):
        tmp = sp_lines[i].split()
        sp[tmp[0]] = i
        sp_num[tmp[0]] = int(tmp[3])

    int_d = {}    # interaction dict
    for i in range(Ni):
        tmp = int_lines[i].split()
        int_d[tmp[0] + " " + tmp[1]] = tmp[3]

    int_mat = np.zeros((Nb, Nb))  # interaction matrix
    for k, v in int_d.items():
        i, j = [sp[l] for l in k.split()]
        int_mat[i, j] = v
        int_mat[j, i] = v

    return int_mat


def parse_control():
    """Return volume"""
    try:
        f = open("CONTROL", "r").readlines()
    except FileNotFoundError:
        sys.exit("CONTROL file not found in the current dir.")
    for line in f:
        if "VOLUME" in line.upper():
            tmp = line.split()
            if len(tmp) == 2:
                L = float(tmp[-1])**(1/3)
                return np.array([L, L, L])
            if len(tmp) == 4:
                return np.array(map(float, tmp[1:4]))


@jit(float64(float64, float64, float64), nopython=True)
def force(nr, a, rc=1.0):
    """
    nr: norm of vector
    a: interaction parameter
    rc: cutoff
    """
    return a * (1 - nr / rc) if nr < rc else 0.0


@jit(float64(float64[:]), nopython=True)
def norm_numba(r):
    rn = 0.0
    for ri in r:
        rn += ri*ri
    return sqrt(rn)


@jit(float64(int64[:], float64[:, :], float64[:, :], float64[:, :], float64))
def virial_p(names, xyz, coeffs, box, rc=1.0):
    N = len(names)
    inv_box = pinv(box)
    V = np.prod(np.diag(box))

    G, Gn = np.zeros(3), np.zeros(3)
    dr, drn = np.zeros(3), np.zeros(3)

    p = 0.0
    for i in range(N):
        for j in range(i):
            dr = xyz[i] - xyz[j]
            G = inv_box @ dr
            Gn = G - np.round(G)
            drn = box @ Gn
            nr = norm_numba(drn)

#            F = force(nr, coeffs[names[i], names[j]], rc)
#            p += nr * F
            pp = coeffs[names[i], names[j]] * (1 - nr / rc) * nr \
                    if nr < rc else 0.0
            p += pp
    return p / (3 * V)


def virial_p_lc(names, xyz, coeffs, box, rc=1.0):
    N = len(names)
    L = np.diag(box)
    V = np.prod(L)
    inv_box = pinv(box)

    Nx = int(max(L // rc))
    lc = LinkCells(L, Nx)
    lc.cell_pairs()
    lc.allocate_atoms(xyz)

    Nintra = sum([len(v) * (len(v) - 1) // 2 for v in lc.cells.values()])
    Ninter = sum([len(lc.cells[p[0]]) * len(lc.cells[p[1]]) for p in lc.pairs])
    dv = np.zeros(Nintra + Ninter)
    G, Gn = np.zeros(3), np.zeros(3)
    dd, ddn = np.zeros(3), np.zeros(3)
    p = 0.0

    for cell in lc.cells.values():
        Na = len(cell)
        for i in range(Na):
            for j in range(i):
                nr = norm_numba(xyz[cell[i]] - xyz[cell[j]])
                F = force(nr, coeffs[names[i], names[j]], rc)
                p += nr * F

    for pair in lc.pairs:
       for i in lc.cells[pair[0]]:
           for j in lc.cells[pair[1]]:
                dd = xyz[i] - xyz[j]
                G = inv_box @ dd
                Gn = G - np.round(G)
                ddn = box @ Gn
                nr = norm_numba(ddn)
                F = force(nr, coeffs[names[i], names[j]], rc)
                p += nr * F

    return p / (3 * V)


if __name__ == "__main__":
    args = docopt(__doc__)
    frame = args["<frame>"]
    A = read_xyzfile(frame)
    names, xyz = A[:, 0], A[:, 1:4]
    names = names.astype(int) - 1
    coeffs = parse_field()
    Ls = parse_control()
    box = np.diag(Ls)
    V = np.prod(Ls)
    N = len(names)
    print("===== Virial pressure calculation for DPD =====")
    print("Box:", Ls)
    print("Interaction coefficients:\n", coeffs)

    print("Calculating pressure...")
    ti = time.time()
    if args["--lc"]:
        pvir = virial_p_lc(names, xyz, coeffs, box)
    else:
        pvir = virial_p(names, xyz, coeffs, box)
    tf = time.time()
    print("Time: %.2f s." % (tf - ti))
    print("Virial pressure:", pvir)
    ptot = pvir + N / V       # kT = 1.0
    print("Total pressure:", ptot)



