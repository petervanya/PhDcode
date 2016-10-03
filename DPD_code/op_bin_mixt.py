#!/usr/bin/env python
"""Usage:
    op_bin_mixt.py normal [--N <N> --rc <rc> --save]
    op_bin_mixt.py lc [--N <N> --rc <rc> --nx <nx>]

[AD HOC] Test order parameter of a binary mixture.
Definition of OP based on Goyal PhD thesis.

Options:
    --N <N>    Number of particles [default: 100]
    --rc <rc>  Cutoff distance [default: 1.3]
    --nx <nx>  Link cells boxes [default: 5]

31/08/16
"""
import numpy as np
from numpy.linalg import norm
from numpy import sqrt
from numba import jit, float64, int64
import itertools
import time
from docopt import docopt
import lmp_lib as ll


def gen_ordered_box(N1, N2):
    types = np.array([1]*N1 + [2]*N2)
    xyz1 = np.random.rand(N1, 3) * [L/2, L, L]
    xyz2 = np.random.rand(N1, 3) * [L/2, L, L]
    xyz2[:, 0] += L/2
    xyz = np.vstack((xyz1, xyz2))
    return types, xyz


def gen_disordered_box(N1, N2):
    types = np.array([1]*N1 + [2]*N2)
    xyz = np.random.rand(N1+N2, 3) * L
    return types, xyz


#def dist_vec(types, xyz, L):
#    N = len(xyz)
#    Np = N * (N-1) // 2
#    dv = np.zeros(Np)                   # dist vec
#    tv = np.zeros((Np, 2), dtype=int)   # type vec
#
#    cell = L * np.eye(3)
#    inv_cell = np.linalg.pinv(cell)
#    dr, drn = np.zeros(3), np.zeros(3)
#    G, Gn = np.zeros(3), np.zeros(3)
#    
#    cnt = 0
#    for i in range(N):
#        for j in range(i):
#            dr = xyz[i] - xyz[j]
#            G = inv_cell @ dr
#            Gn = G - np.round(G)
#            drn = cell @ Gn
#            dv[cnt] = norm_numba(drn)
#            tv[cnt] = [types[i], types[j]]
#            cnt += 1
#    return dv, tv


@jit(float64(float64[:]), nopython=True)
def norm_numba(r):
    rn = 0.0
    for ri in r:
        rn += ri*ri
    return sqrt(rn)


#@jit(float64(int64[:], float64[:, :], float64, float64))
def order_param_naive(types, xyz, L, rc):
    N = len(xyz)
    box = L * np.eye(3)
    inv_box = np.linalg.pinv(box)
    dr, drn = np.zeros(3), np.zeros(3)
    G, Gn = np.zeros(3), np.zeros(3)
    op = []
    
    for i in range(N):
        phi = 0.0
        n1, n2 = 0, 0
        for j in range(N):
            dr = xyz[i] - xyz[j]
            G = inv_box @ dr
            Gn = G - np.round(G)
            drn = box @ Gn
            d = norm(drn)
#            d = norm_numba(drn)
            if i == j:
                continue
            if d < rc:
                if types[j] == 1:
                    n1 += 1
                if types[j] == 2:
                    n2 += 1
        if n1 + n2 != 0:
            phi = (n1 - n2)**2 / (n1 + n2)**2
            op.append(phi)
#        else:
#            print("WARNING: no neighbours within given cutoff.")
#    print(len(op))
    return sum(op) / len(op)


def order_param_lc(types, xyz, lc, L, rc):
    N = len(xyz)
    box = L * np.eye(3)
    inv_box = np.linalg.pinv(box)
    dr, drn = np.zeros(3), np.zeros(3)
    G, Gn = np.zeros(3), np.zeros(3)
    op = []
    Lx = L / Nx
    
    for i in range(N):
        phi = 0.0
        n1, n2 = 0, 0

        id_i = id_from_coord(xyz[i] // Lx, Nx)
        neigh_cells = neighbour_cells(id_i, Nx)
        neigh_atoms = []
        for nc in neigh_cells:
            neigh_atoms.extend(lc[nc])

        for j in neigh_atoms:
            dr = xyz[i] - xyz[j]
            G = inv_box @ dr
            Gn = G - np.round(G)
            drn = box @ Gn
            d = norm(drn)
#            d = norm_numba(drn)
            if i == j:
                continue
            if d < rc:
                if types[j] == 1:
                    n1 += 1
                if types[j] == 2:
                    n2 += 1
        if n1 + n2 != 0:
            phi = (n1 - n2)**2 / (n1 + n2)**2
            op.append(phi)
    return sum(op) / len(op)


# ===== link cells
def cell_coord(id, Nx):
    """Map cell ID to 3d grid number"""
    nx = id // (Nx**2)
    ny = (id - nx * Nx**2) // Nx
    nz = id - nx * Nx**2 - ny * Nx
    return np.array([nx, ny, nz])


def id_from_coord(n, Nx):
    """Map 3d grid number to cell ID"""
    return int(n[0] * Nx**2 + n[1] * Nx + n[2])


def neighbour_cells(id, Nx):
    """Find all neighbouring cells incl the given one"""
    r = cell_coord(id, Nx)
    neighs = []
    tmp = np.arange(3) - 1
    for p in itertools.product(tmp, tmp, tmp):
        neigh = (r + p) % Nx
        neighs.append(neigh)
    return [id_from_coord(neigh, Nx) for neigh in neighs]


def init_link_cells(Lx, Nx):
    """
    Nx: number of cells in one dimension
    Lx: link cell size
    """
    lc = {}
    for i in range(Nx):
        for j in range(Nx):
            for k in range(Nx):
                lc[id_from_coord([i, j, k], Nx)] = []
    return lc


def populate_link_cells(lc, xyz, Lx, Nx):
    """Allocate atoms to link cells. For each atom, find out
    to which cell it belongs and append its name/number to the cell list.
    * lc: link cells"""
    N = len(xyz)
    for i in range(N):
        num = xyz[i] // Lx
        lc[id_from_coord(num, Nx)].append(i)


if __name__ == "__main__":
    args = docopt(__doc__)
    np.random.seed(1234)
    L = 10.0
    f = 0.5
    N = int(args["--N"])
    N1, N2 = int(f * N), int((1-f) * N)
    rc = float(args["--rc"])


    if args["lc"]:
        Nx = int(args["--nx"])
        Lx = L / Nx
        types, xyz = gen_ordered_box(N1, N2)
        if Lx < rc:
            print("WARNING: Link cell size smaller than cutoff radius.")
        print("Testing link cells.")
        lc = init_link_cells(Lx, Nx)
        populate_link_cells(lc, xyz, Lx, Nx)
#        print(lc)
        ti = time.time()
        op = order_param_lc(types, xyz, lc, L, rc)
        print("=== op: %.3f" % op)
        tf = time.time()
        print("Time: %.2f s." % (tf - ti))


    if args["normal"]:
        print("N: %i | L: %.1f | f: %.2f" % (N, L, f))
        # order
        print("Testing ordered box. rc: %.2f" % rc)
        ti = time.time()
        types, xyz = gen_ordered_box(N1, N2)
        op = order_param_naive(types, xyz, L, rc)
        print("=== op: %.3f" % op)
        tf = time.time()
        print("Time: %.2f s." % (tf - ti))
        
        if args["--save"]:
            ll.save_xyzfile("order.xyz", np.vstack((types, xyz.T)).T)
       
        # disorder
        print("Testing disordered box. rc: %.2f" % rc)
        ti = time.time()
        types, xyz = gen_disordered_box(N1, N2)
        op = order_param_naive(types, xyz, L, rc)
        print("=== op: %.3f" % op)
        tf = time.time()
        print("Time: %.2f s." % (tf - ti))
       
        if args["--save"]:
            ll.save_xyzfile("disorder.xyz", np.vstack((types, xyz.T)).T)


