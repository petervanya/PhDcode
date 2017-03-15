#!/usr/bin/env python
"""Usage: surface_tension.py <frames> <coord> [--bins <nb>]

Compute surface tension from Kirkwood relation.

Options:
    --bins <nb>    Number of bins [default: 10]

pv278@cam.ac.uk, 10/03/17
"""
import numpy as np
from numpy.linalg import norm, pinv
from math import sqrt
from numba import jit, float64, int64
from scipy.integrate import simps
import sys, time, glob
from docopt import docopt
from link_cells import LinkCells
from dlms_lib import read_xyzfile2


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
                print(L)
                return np.array([L, L, L])
            if len(tmp) == 4:
                L = np.array(list(map(float, tmp[1:4])))
                print(L)
                return L


@jit(float64(float64, float64, float64), nopython=True)
def force(nr, a, rc=1.0):
    """
    nr: norm of vector
    a: interaction parameter
    rc: cutoff
    """
    return a * (1 - nr / rc) if nr < rc else 0.0


@jit(float64[:](float64[:], float64, float64), nopython=True)
def force_vec(r, a, rc=1.0):
    """
    r: 3d vector
    a: interaction parameter
    rc: cutoff
    """
#    nr = norm_numba(r)
#    return a * (1 - nr / rc) * r / nr if nr < rc else 0.0
    nr = 0.0
    for ri in r:
        nr += ri * ri
    nr = sqrt(nr)
    F = np.zeros(3)
    for i in range(3):
        F[i] = a * (1 - nr / rc) * r[i] / nr if nr < rc else 0.0
    return F


@jit(float64(float64[:]), nopython=True)
def norm_numba(r):
    rn = 0.0
    for ri in r:
        rn += ri * ri
    return sqrt(rn)


@jit(float64(float64), nopython=True)
def theta(x):
    return 1.0 if x > 0.0 else 0.0


@jit(float64(float64, float64, float64), nopython=True)
def ksi(z, zi, zj):
    """Ksi function from Ghoufi, PRE, 2010
    acting as a delta function."""
    zij = zj - zi
    return theta((z - zi) / zij) * theta((zj - z) / zij) / abs(zij)


@jit(float64(float64, float64, float64), nopython=True)
def ksi2(z, zi, zj):
    """Ksi function derived by Mathematica based on Ghoufi, PRE, 2010"""
    return (theta(z - zj) - theta(z - zi)) / (zi - zj)


def align(u, ui, uj, Lu):
    """Consider PBC for the pressure coordinate"""
    if u < Lu / 2 and ui < Lu / 2: uj -= Lu
    if u < Lu / 2 and uj < Lu / 2: ui -= Lu
    if u > Lu / 2 and ui > Lu / 2: uj += Lu
    if u > Lu / 2 and uj > Lu / 2: ui += Lu
    return ui, uj


def pressure_inhomo(names, xyz, coeffs, box, u, Nb, rc=1.0):
    """
    Compute pressure along a coordinate u in [0, 1, 2].
    * Nb: number of bins
    """
    L = np.diag(box)
    V = np.prod(L)
    inv_box = pinv(box)

#    Nx = int(max(L // rc))
#    lc = LinkCells(L, Nx)
#    lc.cell_pairs()
#    lc.allocate_atoms(xyz)
    G, Gn = np.zeros(3), np.zeros(3)
    dr, drn = np.zeros(3), np.zeros(3)

    up = list(set(range(3)).difference([u])) # coords perp. to u
    A = np.prod(L[up])                       # cross-sectional area
    du = L[u] / Nb
    us = np.arange(du / 2, L[u], du)
    P = np.zeros((Nb, 3))

    for ni in range(Nb):
        nums = np.arange(len(xyz))
        nums = nums[(xyz[:, u] > us[ni] - rc / 2) & (xyz[:, u] < us[ni] + rc / 2)]
        Nl = len(nums)
        for i in range(Nl):
            for j in range(i):
                dr = xyz[i] - xyz[j]
                G = inv_box @ dr
                Gn = G - np.round(G)
                drn = box @ Gn
#                if np.any(dr != drn):
#                    ui, uj = align(us[ni], xyz[i, u], xyz[j, u], L[u])
                F = force_vec(drn, coeffs[names[i], names[j]], rc)
                P[ni, :] += drn * F * ksi(us[ni], xyz[i, u], xyz[j, u])
    return P / A



if __name__ == "__main__":
    args = docopt(__doc__)
    frames = sorted(glob.glob(args["<frames>"]))
    Nf = len(frames)
    if Nf == 0:
        sys.exit("No frames captured.")
    coeffs = parse_field()
    Ls = parse_control()
    box = np.diag(Ls)
    V = np.prod(Ls)

    u = int(args["<coord>"])
    Nb = int(args["--bins"])
    up = list(set(range(3)).difference([u])) # coords perp. to u
    if u not in range(3):
        sys.exit("Coordinate must be in [0, 1, 2].")
    du = Ls[u] / Nb
    us = np.arange(du / 2, Ls[u], du)
    P, Pi = np.zeros((Nb, 3)), np.zeros((Nb, 3))
    Pt = np.zeros(Nb)
    gamma, gammai = 0.0, 0.0

    print("===== Surface tension =====")
    print("Pressure along coordinate: %i" % u)
    print("Box: %s" % Ls)
    print("Coeffs:\n", coeffs)

    ti = time.time()
    for frame in frames:
        print("Computing frame %s... " % frame, end="")
        names, xyz = read_xyzfile2(frame)
        names = names.astype(int) - 1
        Pi = pressure_inhomo(names, xyz, coeffs, box, u, Nb)
        gammai = simps(Pi[:, u] - (Pi[:, up[0]] + Pi[:, up[1]]) / 2.0, us)
        gamma += gammai / Nf
        P += Pi / Nf
        print("Gamma:", gammai)
    tf = time.time()
    print("Time: %.2f s." % (tf - ti))

#    print(P)
    Pt = P[:, u] - (P[:, up[0]] + P[:, up[1]]) / 2.0
    np.savetxt("pressure.out", np.c_[us, Pt])
    print("Surface tension:", gamma)


