#!/usr/bin/env python
"""Usage: 
    surf_tension_mdpd.py <frames> <coord> [--bins <nb>]

Compute surface tension using Irving-Kirkwood method.

Arguments:
    <coord>      Coordinate along which to calculate s.t. 'x', 'y', 'z'

Options:
    --bins <nb>  Number of bins [default: 10]

pv278@cam.ac.uk, 13/12/17
"""
import numpy as np
from scipy.integrate import simps
import sys, time, glob
from docopt import docopt
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
            Nbt = int(f[line].split()[-1])
            sp_lines = f[line+1 : line+Nbt+1]
        if "INTERACTIONS" in f[line].upper():
            Ni = int(f[line].split()[-1])
            int_lines = f[line+1 : line+Ni+1]

    sp = {}       # species
    sp_num = {}   # species dict with numbers
    for i in range(Nbt):
        tmp = sp_lines[i].split()
        sp[tmp[0]] = i
        sp_num[tmp[0]] = int(tmp[3])

    int_d = {}    # interaction dict
    for i in range(Ni):
        tmp = int_lines[i].split()
        if tmp[2].lower() != "mdpd":
            sys.exit("Interaction type should be mdpd.")
        int_d[tmp[0] + " " + tmp[1]] = tmp[3:5]

    int_A = np.zeros((Nbt, Nbt))  # interaction matrix
    int_B = np.zeros((Nbt, Nbt))
    for k, v in int_d.items():
        i, j = [sp[l] for l in k.split()]
        int_A[i, j] = v[0]
        int_A[j, i] = v[0]
        int_B[i, j] = v[1]
        int_B[j, i] = v[1]

    return int_A, int_B


def perodic_1D_difference(v1, v2, Lu):
    outer_prod = np.subtract.outer(v1, v2)
    outer_prod = np.where((np.abs(outer_prod) > Lu / 2), \
            (Lu - np.abs(outer_prod)) * -1 * np.sign(outer_prod), \
            outer_prod)
    return outer_prod


def distance_matrix(xyz1, xyz2, L):
    """need to get in the periodic boundary"""
    x = perodic_1D_difference(xyz1[:, 0], xyz2[:, 0], L[0])
    y = perodic_1D_difference(xyz1[:, 1], xyz2[:, 1], L[1])
    z = perodic_1D_difference(xyz1[:, 2], xyz2[:, 2], L[2])
    return x, y, z


def local_density(r, rd, rc=1.0):
    """r: (N, N) istance matrix after PBC"""
    rho = 15.0 / (2.0 * np.pi * rd**3) * (r < rd) * (r > 0) * (1 - r / rd)**2
    return np.sum(rho, 1)


def boxsplit(xyz, n, ui, u, L, rc=1.0):
    """ n: list of bead numbers"""
    if (ui >= rc) and (ui <= (L[u] - rc)):
        xyz1 = xyz[np.logical_and((xyz[:, u] > ui - rc), (xyz[:, u] < ui))]
        xyz2 = xyz[np.logical_and((xyz[:, u] < ui + rc), (xyz[:, u] > ui))]
        n1 = n[np.logical_and((xyz[:, u] > ui - rc), (xyz[:, u] < ui))]
        n2 = n[np.logical_and((xyz[:, u] < ui + rc), (xyz[:, u] > ui))]
    elif ui < rc:
        xyz1 = xyz[np.logical_or((xyz[:, u] > L[u] - rc), (xyz[:, u] < ui))]
        xyz2 = xyz[np.logical_and((xyz[:, u] < ui + rc), (xyz[:, u] > ui))]
        n1 = n[np.logical_or((xyz[:, u] > L[u] - rc), (xyz[:, u] < ui))]
        n2 = n[np.logical_and((xyz[:, u] < ui + rc), (xyz[:, u] > ui))]
    elif ui > L[u] - rc:
        xyz1 = xyz[np.logical_and((xyz[:, u] > ui - rc), (xyz[:, u] < ui))]
        xyz2 = xyz[np.logical_or((xyz[:, u] < rc), (xyz[:, u] > ui))]
        n1 = n[np.logical_and((xyz[:, u] > ui - rc), (xyz[:, u] < ui))]
        n2 = n[np.logical_or((xyz[:, u] < rc), (xyz[:, u] > ui))]
    return n1, xyz1, n2, xyz2


def pressure_tensor(nm, xyz, rho, cA, cB, Ls, u, Nb, rd=0.75, rc=1.0):
    us = np.linspace(0, Ls[u], Nb+1)[:-1]
    P = np.zeros((Nb, 3))
    n = np.arange(len(xyz))

    for nb in range(Nb):
        nl, xyzl, nr, xyzr = boxsplit(xyz, n, us[nb], u, Ls)
        dms = distance_matrix(xyzl, xyzr, Ls)
        r = np.sqrt(dms[0]**2 + dms[1]**2 + dms[2]**2)
        
        Amat = np.array([[cA[nm[i], nm[j]] for i in nl] for j in nr]).T
        Bmat = np.array([[cB[nm[i], nm[j]] for i in nl] for j in nr]).T

        F = Amat * (1 - r) * (r < 1) * (r > 0) + \
                Bmat * (1 - r / rd) * np.add.outer(rho[nl], rho[nr]) * \
                (r < rd) * (r > 0)

        Fx = F * dms[0] / r
        Fy = F * dms[1] / r
        Fz = F * dms[2] / r
	
        p_xx = np.sum(dms[0] / np.abs(dms[u]) * Fx)
        p_yy = np.sum(dms[1] / np.abs(dms[u]) * Fy)
        p_zz = np.sum(dms[2] / np.abs(dms[u]) * Fz)
        P[nb] = [p_xx, p_yy, p_zz]
    print("Pressure:\n",P)
    return P


if __name__ == "__main__":
    args = docopt(__doc__)
    frames = sorted(glob.glob(args["<frames>"]))
    Nf = len(frames)
    if Nf == 0:
        sys.exit("No frames captured.")
    Ls, rd = parse_control()
    coeffsA, coeffsB = parse_field()

    coord = args["<coord>"]
    coords = {"x": 0, "y": 1, "z": 2}
    if coord not in coords.keys():
        sys.exit("Coordinate must be in %s." % list(coords.keys()))
    u = coords[coord]
    Nb = int(args["--bins"])
    us = np.linspace(0, Ls[u], Nb+1)[:-1]
    uperp = list(set(range(3)).difference([u])) # coords perp. to u
    Area = np.prod(Ls[uperp])                   # cross-sectional area

    P = np.zeros((Nb, 3))
    Pt = np.zeros(Nb)
    gammas = np.zeros(Nf)

    print("===== Surface tension for MDPD =====")
    print("Frames: %i | Box: %s | rd: %.2f" % (Nf, Ls, rd))
    print("Nu: %i | du: %.4f" % (len(us), Ls[u] / len(us)))
    print("Coeffs:\n", coeffsA, "\n", coeffsB)

    ti = time.time()
    for j in range(Nf):
        nm, xyz = read_xyzfile2(frame)
        nm = nm.astype(int) - 1
        dms = distance_matrix(xyz, xyz, Ls)
        r = np.sqrt(dms[0]**2 + dms[1]**2 + dms[2]**2)
        rho = local_density(r, rd)

        Pi = pressure_tensor(nm, xyz, rho, coeffsA, coeffsB, \
                Ls, u, Nb, rd, rc=1.0)
        Pti = (Pi[:, u] - (Pi[:, uperp[0]] + Pi[:, uperp[1]]) / 2.0) / Area
        Pt += Pti / Nf
        gammas[j] = simps(Pti, us) / 2.0
        print(frames[j], " | Gamma = ", gammas[j])

    print("Time: %.2f s." % (time.time() - ti))
    np.savetxt("pressure.out", np.c_[us, Pt])
    gamma_avg = simps(Pt, us) / 2.0
    print("Final gamma: %f | Std: %f" % (gamma_avg, np.std(gammas)))


