#!/usr/bin/env python
"""Usage:
    widom_insertion.py <xyz> <type> <Nins> [--seed <s>]

Calculate the excess chemical potential of a DPD fluid
via Widom particle insertion. No link cells, perhaps add later.

Arguments:
    <xyz>          xyz frame
    <type>         Particle type (number) from those in the frame
    <Nins>         Number of trial insertions

Options:
    --seed <s>     Random seed [default: 49]

19/02/18
"""
import numpy as np
from numpy.linalg import norm
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
        if "volume" in line.lower():
            tmp = line.split()
            if len(tmp) == 2:
                L = float(tmp[-1])**(1/3)
                L = np.ones(3) * L
            elif len(tmp) == 4:
                L = np.array(list(map(float, tmp[1:4])))
    return L


def parse_field():
    """Extract coefficients from the field file. Assume rc = 1.0"""
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
        if tmp[2].lower() != "dpd":
            sys.exit("Interaction type should be dpd.")
        ip_d[tmp[0] + " " + tmp[1]] = tmp[3:5]

    ip_A = np.zeros((Nbt, Nbt))  # interaction matrix
    for k, v in ip_d.items():
        i, j = [sp[l] for l in k.split()]
        ip_A[i, j] = v[0]
        ip_A[j, i] = v[0]

    return ip_A


@jit(float64(float64[:]), nopython=True)
def norm_numba(r):
    rn = 0.0
    for ri in r:
        rn += ri * ri
    return np.sqrt(rn)


@jit(float64(float64), nopython=True)
def wr(nr):
    """Weight function, r>0 is mutual distance"""
    return (1 - nr) if nr < 1.0 else 0.0


@jit(nopython=True)
def excess_pe(xyz, nm, ip, box, inv_box, pos):
    """Excess potential energy of N+1-st particle"""
    N = len(xyz)
    rij, g = np.zeros(3), np.zeros(3)
    U = 0.0

    for i in range(N):
        rij = xyz[i] - pos
        g = inv_box @ rij
        g = g - np.round_(g, 0, np.empty_like(g))
        rij = box @ g
        A = ip[typ-1, nm[i]-1]
        U += A * wr(norm_numba(rij))**2 / 2.0
    return U


@jit
def diff_pe(xyz, nm, ip, box, inv_box, pos):
    """Difference of potential energies of N+1 and N particles.
    Equivalent to excess pot e, but two orders of magnitude slower."""
    xyz1 = np.vstack((xyz, pos))
    nm1 = np.hstack((nm, typ))
    rij, g = np.zeros(3), np.zeros(3)
    U0 = 0.0
    U1 = 0.0

    for i in range(len(xyz)):
        for j in range(i):
            rij = xyz[i] - xyz[j]
            g = inv_box @ rij
            g = g - np.round_(g, 0, np.empty_like(g))
            rij = box @ g
            A = ip[nm[i]-1, nm[j]-1]
            U0 += A * wr(norm_numba(rij))**2 / 2.0

    for i in range(len(xyz1)):
        for j in range(i):
            rij = xyz1[i] - xyz1[j]
            g = inv_box @ rij
            g = g - np.round_(g, 0, np.empty_like(g))
            rij = box @ g
            A = ip[nm1[i]-1, nm1[j]-1]
            U1 += A * wr(norm_numba(rij))**2 / 2.0
    return U1 - U0


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
    ip = parse_field()
    L = parse_control()
    box = np.diag(L)
    inv_box = np.diag(1.0 / L)
    Uex = np.zeros(Nins)

    print("===== Widom particle insertion =====")
    print("Box: %s" % L)
    print("Measuring excess potential energy...")
    ti = time.time()
    for i in range(Nins):
        pos = np.random.rand(3) * L
        Uex[i] = excess_pe(xyz, nm, ip, box, inv_box, pos)
#        Uex[i] = diff_pe(xyz, nm, ip, box, inv_box, pos)
    print("Done. Time: %.2f s" % (time.time() - ti))

    print(Uex)
    muex = -np.log(np.average(np.exp(-Uex)))
    print("Excess chemical potential: %.6f" % muex)
    

