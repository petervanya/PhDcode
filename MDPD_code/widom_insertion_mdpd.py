#!/usr/bin/env python
"""Usage:
    widom_insertion_mdpd.py <xyz> <type> <Nins> [--seed <s>]

Calculate excess chemical potential via Widom particle insertion
for MDPD fluid.
No link cells, add later.

Arguments:
    <xyz>          xyz frame
    <type>         Particle type (number) from those in the frame
    <Nins>         Number of trial insertions

Options:
    --seed <s>     Random seed [default: 49]

20/02/18
"""
import numpy as np
from numpy import pi
from numba import jit, float64
import sys, time
from dlms_lib import read_xyzfile2
from docopt import docopt


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


@jit(float64(float64[:]), nopython=True)
def norm_numba(r):
    rn = 0.0
    for ri in r:
        rn += ri * ri
    return np.sqrt(rn)


@jit(float64(float64, float64), nopython=True)
def wr(nr, rd):
    """Weight function, r>0 is mutual distance"""
    return (1 - nr/rd) if nr < rd else 0.0


@jit(nopython=True)
def excess_pe(xyz, nm, ip_A, ip_B, rd, box, inv_box, pos):
    """Excess potential energy of N+1-st particle"""
    N = len(xyz)
    rij, g = np.zeros(3), np.zeros(3)
    rhoA, rhoB = 0.0, 0.0

    for i in range(N):
        rij = xyz[i] - pos
        g = inv_box @ rij
        g = g - np.round_(g, 0, np.empty_like(g))
        rij = box @ g
        nr = norm_numba(rij)
        A = ip_A[typ-1, nm[i]-1]
        B = ip_B[typ-1, nm[i]-1]
        rhoA += 15/(2*pi) * wr(nr, 1.0)**2
        rhoB += 15/(2*pi*rd**3) * wr(nr, rd)**2
    U = 2 * pi/30 * A * rhoA + pi*rd**4/30 * B * rhoB**2
    return U


@jit(nopython=True)
def local_density(xyz, rd, box, inv_box):
    N = len(xyz)
    rij, g = np.zeros(3), np.zeros(3)
    drho = 0.0
    rho = np.zeros(N)

    for i in range(N):
        for j in range(i):
            rij = xyz[i] - xyz[j]
            g = inv_box @ rij
            g = g - np.round_(g, 0, np.empty_like(g))
            rij = box @ g
            nr = norm_numba(rij)
            drho = 15/(2*pi*rd**3) * wr(nr, rd)**2
            rho[i] += drho
            rho[j] += drho
    return rho


if __name__ == "__main__":
    args = docopt(__doc__)
    seed = int(args["--seed"])
    np.random.seed(seed)
    frame = args["<xyz>"]
    nm, xyz = read_xyzfile2(frame)
    typs = set(nm)
    typ = int(args["<type>"])
    if typ not in typs:
        sys.exit("Requested particle type not found in the frame.")
    Nins = int(args["<Nins>"])
    ip_A, ip_B = parse_field()
    L, rd = parse_control()
    box = np.diag(L)
    inv_box = np.diag(1.0 / L)
    Uex = np.zeros(Nins)

    print("===== Widom particle insertion =====")
    print("Box: %s" % L)
    print("Manybody cutoff:", rd)
    print("Parameters:\n", ip_A, "\n", ip_B)
    print("Measuring excess potential energy...")

    rhoA = local_density(xyz, 1.0, box, inv_box)
    rhoB = local_density(xyz, rd, box, inv_box)

    ti = time.time()
    for i in range(Nins):
        pos = np.random.rand(3) * L
#        Uex[i] = excess_pe(xyz, nm, ip_A, ip_B, rd, box, inv_box, pos)

        xyz1 = np.vstack((xyz, pos))
        rhoA1 = local_density(xyz1, 1.0, box, inv_box)
        rhoB1 = local_density(xyz1, rd, box, inv_box)
 
        A, B = ip_A[0, 0], ip_B[0, 0]
        U0 = pi/30 * A * np.sum(rhoA) + pi*rd**4/30 * B * np.sum(rhoB**2)
        U1 = pi/30 * A * np.sum(rhoA1) + pi*rd**4/30 * B * np.sum(rhoB1**2)
        Uex[i] = U1 - U0
    print("Done. Time: %.2f s" % (time.time() - ti))

    print(Uex)
    muex = -np.log(np.average(np.exp(-Uex)))
    print("Excess chemical potential: %.6f" % muex)
    

