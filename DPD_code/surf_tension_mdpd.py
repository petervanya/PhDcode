#!/usr/bin/env python
"""Usage: 
    surf_tension_mdpd.py <frames> <coord> (kb | ik) [--bins <nb>]

Compute surface tension from Kirkwood relation.

Arguments:
    <coord>      Coordinate along which to calculate s.t. 'x', 'y', 'z'
    kb           Kirkwood-Buff version
    ik           Irving-Kirkwood version, with pressure along coordinate

Options:
    --bins <nb>  Number of bins [default: 10]

pv278@cam.ac.uk, 22/03/17
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


def parse_control():
    """Return volume"""
    rd = 0.0
    try:
        f = open("CONTROL", "r").readlines()
    except FileNotFoundError:
        sys.exit("CONTROL file not found in the current dir.")
    for line in f:
        if "manybody cutoff" in line.lower():
            rd = float(line.split()[-1])
        if "volume" in line.lower():
            tmp = line.split()
            if len(tmp) == 2:
                L = float(tmp[-1])**(1/3)
                L = np.array([L, L, L])
            elif len(tmp) == 4:
                L = np.array(list(map(float, tmp[1:4])))
    return L, rd


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
        if tmp[2] != "mdpd":
            sys.exit("Interaction type: %s. Should be: mdpd." % tmp[1])
        int_d[tmp[0] + " " + tmp[1]] = tmp[3:5]

    int_A = np.zeros((Nb, Nb))  # interaction matrix
    int_B = np.zeros((Nb, Nb))
    for k, v in int_d.items():
        i, j = [sp[l] for l in k.split()]
        int_A[i, j] = v[0]
        int_A[j, i] = v[0]
        int_B[i, j] = v[1]
        int_B[j, i] = v[1]

    return int_A, int_B


@jit(float64(float64, float64), nopython=True)
def weight_rho(nr, rd=0.75):
    return 15.0 / (2 * np.pi * rd**3) * (1.0 - nr / rd)**2 if nr < rd else 0.0


def density(names, xyz, Ls, rd=0.75, rc=1.0):
    """Now works with only one bead type!"""
    N = len(xyz)
    box = np.diag(Ls)  # FIX
    inv_box = pinv(box)

    Nx = int(max(Ls // rc))
    lc = LinkCells(Ls, Nx)
    lc.cell_pairs()
    lc.allocate_atoms(xyz)

    G, Gn = np.zeros(3), np.zeros(3)
    dr, drn = np.zeros(3), np.zeros(3)
    rho = np.zeros(N)

    for cell in lc.cells.values():
        Na = len(cell)
        for i in range(Na):
            for j in range(i):
                nr = norm_numba(xyz[cell[i]] - xyz[cell[j]])
                rho[i] += weight_rho(nr, rd)
                rho[j] += weight_rho(nr, rd)

    for pair in lc.pairs:
       for i in lc.cells[pair[0]]:
           for j in lc.cells[pair[1]]:
                dd = xyz[i] - xyz[j]
                G = inv_box @ dd
                Gn = G - np.round(G)
                ddn = box @ Gn
                nr = norm_numba(ddn)
                rho[i] += weight_rho(nr, rd)
                rho[j] += weight_rho(nr, rd)

    return rho


@jit(float64(float64, float64), nopython=True)
def weight(nr, rc=1.0):
    return 1.0 - nr / rc if nr < rc else 0.0


@jit(float64[:](float64[:], float64, float64, float64, float64,\
        float64, float64), nopython=True)
def force_vec_mdpd(r, A, B, rho_i, rho_j, rd=0.75, rc=1.0):
    nr = 0.0
    for ri in r:
        nr += ri * ri
    nr = sqrt(nr)
    F = np.zeros(3)
    for i in range(3):
        F[i] = (A * weight(nr, rc) + B * (rho_i + rho_j) * weight(nr, rd)) \
                * r[i] / nr
    return F


@jit(float64(float64, float64, float64, float64, float64, \
        float64, float64), nopython=True)
def force_mdpd(nr, A, B, rho_i, rho_j, rd=0.75, rc=1.0):
    """
    nr: norm of vector
    A, B: interaction parameters
    rho_i, rho_j: local densities
    rc: cutoff
    """
    return A * weight(nr, rc) + B * (rho_i + rho_j) * weight(nr, rd)


@jit(float64[:](float64[:], float64, float64), nopython=True)
def force_vec(r, a, rc=1.0):
    """
    r: 3d vector
    a: interaction parameter
    rc: cutoff
    """
    nr = 0.0
    for ri in r:
        nr += ri * ri
    nr = sqrt(nr)
    F = np.zeros(3)
    for i in range(3):
        F[i] = a * (1 - nr / rc) * r[i] / nr if nr < rc else 0.0
    return F


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


def add_images(xyz, Ls, u, rc=1.0):
    """Transfer images from one end of the box to the other,
    along a coord u in [0, 1, 2]"""
    xyz_low = xyz[xyz[:, u] < rc]
    xyz_low[:, u] += Ls[u]
    xyz_high = xyz[xyz[:, u] > Ls[u] - rc]
    xyz_high[:, u] += Ls[u]
    return np.r_[xyz, xyz_low, xyz_high]


def pressure_vec(names, xyz, rho, cA, cB, box, u, Nbins, rd=0.75, rc=1.0):
    """
    Compute pressure along a coordinate u in [0, 1, 2].
    * Nbins: number of bins
    """
    Ls = np.diag(box)
    V = np.prod(Ls)
    inv_box = pinv(box)

    G, Gn = np.zeros(3), np.zeros(3)
    dr, drn = np.zeros(3), np.zeros(3)

    uperp = list(set(range(3)).difference([u])) # coords perp. to u
    Area = np.prod(Ls[uperp])                   # cross-sectional area
    du = Ls[u] / Nbins
    us = np.arange(du / 2, Ls[u], du)
    P = np.zeros((Nbins, 3))

    for ni in range(Nbins):
        nums = np.arange(len(xyz))
        nums = nums[(xyz[:, u] > us[ni] - rc / 2) & (xyz[:, u] < us[ni] + rc / 2)]
        Nl = len(nums)
        for i in range(Nl):
            for j in range(i):
                dr = xyz[nums[i]] - xyz[nums[j]]
                G = inv_box @ dr
                Gn = G - np.round(G)
                drn = box @ Gn
                drn[u] = dr[u]

                A = cA[names[i], names[j]]
                B = cB[names[i], names[j]]
                F = force_vec_mdpd(drn, A, B, rho[i], rho[j], rd, rc)
                P[ni, :] += drn * F * ksi(us[ni], xyz[i, u], xyz[j, u])
    return P / Area


def surf_tension_kb(names, xyz, rho, cA, cB, box, u, Nbins, rd=0.75, rc=1.0):
    """Kirkwood-Buff version of computing surface tension"""
    N = len(xyz)
    L = np.diag(box)
    V = np.prod(L)
    inv_box = pinv(box)

    G, Gn = np.zeros(3), np.zeros(3)
    dr, drn = np.zeros(3), np.zeros(3)

    uperp = list(set(range(3)).difference([u])) # coords perp. to u
    Area = np.prod(L[uperp])                    # cross-sectional area
    gamma = 0.0
    for i in range(N):
        for j in range(i):
            dr = xyz[i] - xyz[j]
            G = inv_box @ dr
            Gn = G - np.round(G)
            drn = box @ Gn

            A = cA[names[i], names[j]]
            B = cB[names[i], names[j]]
            nr = norm_numba(drn)
            F = force_mdpd(nr, A, B, rho[i], rho[j], rd, rc)
            gamma += (3 * drn[u]**2 - nr**2) / (2 * nr) * F
    return gamma / (2 * Area)


def surf_tension_kb_lc(names, xyz, rho, cA, cB, box, u, Nbins, rd=0.75, rc=1.0):
    """Kirkwood-Buff version of computing surface tension"""
    N = len(xyz)
    Ls = np.diag(box)
    V = np.prod(Ls)
    inv_box = pinv(box)

    uperp = list(set(range(3)).difference([u])) # coords perp. to u
    Area = np.prod(Ls[uperp])                   # cross-sectional area

    Nx = int(max(Ls // rc))
    lc = LinkCells(Ls, Nx)
    lc.cell_pairs()
    lc.allocate_atoms(xyz)

    G, Gn = np.zeros(3), np.zeros(3)
    dr, drn = np.zeros(3), np.zeros(3)

    gamma = 0.0
    for cell in lc.cells.values():
        Na = len(cell)
        for i in range(Na):
            for j in range(i):
                drn = xyz[cell[i]] - xyz[cell[j]]
                A = cA[names[cell[i]], names[cell[j]]]
                B = cB[names[cell[i]], names[cell[j]]]
                nr = norm_numba(drn)
                F = force_mdpd(nr, A, B, rho[cell[i]], rho[cell[j]], rd, rc)
                gamma += (3 * drn[u]**2 - nr**2) / (2 * nr) * F

    for pair in lc.pairs:
       for i in lc.cells[pair[0]]:
           for j in lc.cells[pair[1]]:
                dr = xyz[i] - xyz[j]
                G = inv_box @ dr
                Gn = G - np.round(G)
                drn = box @ Gn
 
                A = cA[names[i], names[j]]
                B = cB[names[i], names[j]]
                nr = norm_numba(drn)
                F = force_mdpd(nr, A, B, rho[i], rho[j], rd, rc)
                gamma += (3 * drn[u]**2 - nr**2) / (2 * nr) * F
    return gamma / (2 * Area)


if __name__ == "__main__":
    args = docopt(__doc__)
    frames = sorted(glob.glob(args["<frames>"]))
    Nf = len(frames)
    if Nf == 0:
        sys.exit("No frames captured.")
    coeffsA, coeffsB = parse_field()
    Ls, rd = parse_control()
    box = np.diag(Ls)
    V = np.prod(Ls)

    coord = args["<coord>"]
    coords = {"x": 0, "y": 1, "z": 2}
    if coord not in coords.keys():
        sys.exit("Coordinate must be in ['x', 'y', 'z'].")
    u = coords[coord]
    Nbins = int(args["--bins"])
    up = list(set(range(3)).difference([u])) # coords perp. to u
    du = Ls[u] / Nbins
    us = np.arange(du / 2, Ls[u], du)
    P, Pi = np.zeros((Nbins, 3)), np.zeros((Nbins, 3))
    Pt = np.zeros(Nbins)
    gamma, gammai = 0.0, 0.0

    print("===== Surface tension =====")
    print("Pressure along coordinate: %s | method: %s" % \
            (coord, "IK" if args["ik"] else "KB"))
    print("Frames: %i | Box: %s | rd: %.2f" % (Nf, Ls, rd))
    print("Coeffs:\n", coeffsA, "\n", coeffsB)

    if args["ik"]:
        ti = time.time()
        for frame in frames:
            print("Computing frame %s... " % frame, end="")
            names, xyz = read_xyzfile2(frame)
            names = names.astype(int) - 1
            rho = density(names, xyz, Ls, rd, rc=1.0)
            xyz = add_images(xyz, Ls, u)
            Pi = pressure_vec(names, xyz, rho, coeffsA, coeffsB, \
                    box, u, Nbins, rd, rc=1.0)
            gammai = simps(Pi[:, u] - (Pi[:, up[0]] + Pi[:, up[1]]) / 2.0, us)
            gamma += gammai / Nf
            P += Pi / Nf
            print("Gamma:", gammai)
        tf = time.time()
        print("Time: %.2f s." % (tf - ti))
      
        Pt = P[:, u] - (P[:, up[0]] + P[:, up[1]]) / 2.0
        np.savetxt("pressure.out", np.c_[us, Pt])
        print("Surface tension:", gamma)

    if args["kb"]:
        ti = time.time()
        for frame in frames:
            print("Computing frame %s... " % frame, end="")
            names, xyz = read_xyzfile2(frame)
            names = names.astype(int) - 1
            rho = density(names, xyz, Ls, rd, rc=1.0)
            gammai = surf_tension_kb(names, xyz, rho, coeffsA, coeffsB, \
                    box, u, Nbins, rd, rc=1.0)
            gamma += gammai / Nf
            print("Gamma:", gammai)
        tf = time.time()
        print("Time: %.2f s." % (tf - ti))
      
        print("Surface tension:", gamma)


