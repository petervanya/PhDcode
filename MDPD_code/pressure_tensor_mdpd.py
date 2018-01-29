#!/usr/bin/env python
"""Usage:
    pressure_tensor_mdpd.py <frames> <coord>

Compute the pressure tensor and the surface tension
for a several xyz frames.

pv278@cam.ac.uk, 15/01/18
"""
import numpy as np
from numba import jit
from scipy.integrate import simps
import sys, time, glob
from docopt import docopt
from dlms_lib import read_xyzfile2


def sort_frame_by(frame):
    return int(frame.split("/")[-1].rstrip(".xyz").split("_")[-1])


def parse_control():
    """Return volume"""
    rd = 1.0
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
            ip_lines = f[line+1 : line+Ni+1]

    sp = {}       # species
    sp_num = {}   # species dict with numbers
    for i in range(Nbt):
        tmp = sp_lines[i].split()
        sp[tmp[0]] = i
        sp_num[tmp[0]] = int(tmp[3])

    ip_d = {}    # interaction dict
    for i in range(Ni):
        tmp = ip_lines[i].split()
        if tmp[2].lower() != "mdpd":
            sys.exit("Interaction type should be mdpd.")
        ip_d[tmp[0] + " " + tmp[1]] = tmp[3:5]

    ip_A = np.zeros((Nbt, Nbt))  # interaction matrix
    ip_B = np.zeros((Nbt, Nbt))
    for k, v in ip_d.items():
        i, j = [sp[l] for l in k.split()]
        ip_A[i, j] = v[0]
        ip_A[j, i] = v[0]
        ip_B[i, j] = v[1]
        ip_B[j, i] = v[1]

    return ip_A, ip_B


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


def pressure_tensor(nm, xyz, rho, cA, cB, Ls, rd=0.75):
    n = np.arange(len(xyz))
    V = np.prod(Ls)

    Dx, Dy, Dz = distance_matrix(xyz, xyz, Ls)
    r = np.sqrt(Dx[0]**2 + Dy[1]**2 + Dz[2]**2)

    Amat = np.array([[cA[nm[i], nm[j]] for i in n] for j in n]).T
    Bmat = np.array([[cB[nm[i], nm[j]] for i in n] for j in n]).T

    Fnorm = Amat * (1 - r) * (r < 1) * (r > 0) + \
            Bmat * (1 - r / rd) * np.add.outer(rho, rho) * \
            (r < rd) * (r > 0)

    Fx = Fnorm * Dx[0] / r
    Fy = Fnorm * Dx[1] / r
    Fz = Fnorm * Dx[2] / r

    return np.array([np.sum(Fx*Dx), np.sum(Fy*Dy), np.sum(Fz*Dz)]) / (3*V)


@jit(nopython=True)
def wr(nr, rd):
    if nr < rd:
        return 1.0 - nr
    else:
        return 0

@jit(nopython=True)
def pressure_tensor2(nm, xyz, rho, cA, cB, Ls, rd=0.75):
    P = np.zeros(3)
    V = np.prod(Ls)
    N = len(xyz)
    box = np.diag(Ls)
    inv_box = np.linalg.pinv(box)
    rij, G = np.zeros(3), np.zeros(3)

    for i in range(N):
        for j in range(i):
            rij = xyz[i] - xyz[j]
            G = inv_box @ rij
            G = G - np.round_(G, 0, np.empty_like(G))
            rij = box @ G
            nr = np.linalg.norm(rij)
            Fij = rij / nr * (cA[nm[i], nm[j]] * wr(nr, 1.0) + \
                    cB[nm[i], nm[j]] * wr(nr, rd) * (rho[i] + rho[j]))
            P += Fij * rij / (3 * V)
    return P


if __name__ == "__main__":
    args = docopt(__doc__)
    frames = glob.glob(args["<frames>"])
    frames.sort(key=sort_frame_by)
    Nf = len(frames)
    if Nf == 0:
        sys.exit("No frames captured.")
    Ls, rd = parse_control()
    ip_A, ip_B = parse_field()

    coord = args["<coord>"]
    coords = {"x": 0, "y": 1, "z": 2}
    if coord not in coords.keys():
        sys.exit("Coordinate must be in %s." % list(coords.keys()))
    u = coords[coord]
    uperp = list(set(range(3)).difference([u])) # coords perp. to u
    Area = np.prod(Ls[uperp])

    P = np.zeros((Nf, 3))
    gammas = np.zeros(Nf)

    print("===== Pressure tensor for MDPD =====")
    print("Frames: %i | Box: %s | rd: %.2f" % (Nf, Ls, rd))
    print("Coeffs:\n", ip_A, "\n", ip_B)
    print("Box length: %.1f | Cross-area: %.1f" % (Ls[u], Area))

    ti = time.time()
    for j in range(Nf):
        nm, xyz = read_xyzfile2(frames[j])
        nm -= 1
        dms = distance_matrix(xyz, xyz, Ls)
        r = np.sqrt(dms[0]**2 + dms[1]**2 + dms[2]**2)
        rho = local_density(r, rd)

        P[j] = pressure_tensor2(nm, xyz, rho, ip_A, ip_B, Ls, rd)
        gammas[j] = (P[j, u] - (P[j, uperp[0]] + P[j, uperp[1]]) / 2.0) * Ls[u] / 2
        print("%s | Pvir: %7.3f | Pdiag: %s | Gamma: %s " % \
                (frames[j], np.sum(P[j]), P[j], gammas[j]))

    print("Time: %.2f s." % (time.time() - ti))
    print("Final gamma: %f | Std: %f" % (np.average(gammas), np.std(gammas)))
