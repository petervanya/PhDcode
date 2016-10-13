#!/usr/bin/env python
"""Usage:
    op_bin_mixt.py <frames> [--rc <rc>]
    op_bin_mixt.py test normal [--N <N> --f <f> --rc <rc> --save]
    op_bin_mixt.py test linkcells [--N <N> --f <f> --rc <rc> --nx <nx>]

Compute order parameter of a binary mixture, using link cells.
Definition of OP based on Goyal PhD thesis.

Options:
    --N <N>    Number of particles [default: 100]
    --f <f>    Fraction of A particles [default: 0.5]
    --rc <rc>  Cutoff distance [default: 1.3]
    --nx <nx>  Link cells boxes in one dim [default: 5]

03/10/16 (revised)
"""
import numpy as np
from numpy.linalg import norm
from numpy import sqrt
from numba import jit, float64, int64
import itertools
import time
import glob, sys
from docopt import docopt
from lmp_lib import read_xyzfile


def gen_ordered_box(N1, N2):
    types = np.array([1]*N1 + [2]*N2)
    xyz1 = np.random.rand(N1, 3) * [L/2, L, L]
    xyz2 = np.random.rand(N2, 3) * [L/2, L, L]
    xyz2[:, 0] += L/2
    xyz = np.vstack((xyz1, xyz2))
    return types, xyz


def gen_disordered_box(N1, N2):
    types = np.array([1]*N1 + [2]*N2)
    xyz = np.random.rand(N1+N2, 3) * L
    return types, xyz


def guess_box_size(xyz):
    """Infer box size from xyz matrix"""
    return np.round(np.max(xyz[:, 1] - np.min(xyz[:, 1]), 0))


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
#            d = norm(drn)
            d = norm_numba(drn)
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


def order_param_lc(types, xyz, lc, L, rc):
    N = len(xyz)
    N1, N2 = sum(types == 1), sum(types == 2)
    f1, f2 = N1 / N, N2 / N
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
#            d = norm(drn)
            d = norm_numba(drn)
            if i == j:
                continue
            if d < rc:
                if types[j] == 1:
                    n1 += 1
                if types[j] == 2:
                    n2 += 1
        if n1 + n2 != 0:
            phi = (n1 - n2)**2 / (n1 + n2)**2
#            phi = (n1 - n2 * f1/f2)**2 / (n1 + n2 * f1/f2)**2
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
        num = xyz[i] // Lx % Nx
        lc[id_from_coord(num, Nx)].append(i)


if __name__ == "__main__":
    args = docopt(__doc__)
    rc = float(args["--rc"])

    if not args["test"]:
        print("===== Order parameter =====")
        frames = glob.glob(args["<frames>"])
        Nf = len(frames)
        if Nf == 0:
            sys.exit("No files captured.")
        A = read_xyzfile(frames[0])
        N = len(A)
        L = guess_box_size(A[:, 1:])
        Nx = int(L // rc)    # choose largest possible cell number
        print("Frames: %i | N: %i | Nx: %i | Box: %.2f" % (Nf, N, Nx, L))

        OP = 0.0
        ti = time.time()
        for frame in frames:
            A = read_xyzfile(frame)
            types, xyz = A[:, 0], A[:, 1:]
            L = guess_box_size(xyz)
            Lx = L / Nx
            lc = init_link_cells(Lx, Nx)
            populate_link_cells(lc, xyz, Lx, Nx)
            op = order_param_lc(types, xyz, lc, L, rc)
            print("%.3f   " % op, end="", flush=True)
            OP += op / Nf
        tf = time.time()
        print("\nOrder parameter: %.3f\nTime: %.2f s." % (OP, tf - ti))


    else:
        np.random.seed(1234)
        L = 10.0
        f = float(args["--f"])
        N = int(args["--N"])
        N1, N2 = int(f * N), int((1-f) * N)
        print("====== Testing order parameter script =====")
        print("N: %i | L: %.1f | f: %.2f | rc: %.2f" % (N, L, f, rc))

        if args["linkcells"]:
            Nx = int(args["--nx"])
            Lx = L / Nx
            if Lx < rc:
                print("WARNING: Link cell size smaller than cutoff radius.")
 
            print("Using link cells.\nOrdered box")
            types, xyz = gen_ordered_box(N1, N2)
            lc = init_link_cells(Lx, Nx)
            populate_link_cells(lc, xyz, Lx, Nx)
            ti = time.time()
            op = order_param_lc(types, xyz, lc, L, rc)
            tf = time.time()
            print("=== op: %.3f\nTime: %.2f s." % (op, tf - ti))
 
            print("Disordered box")
            types, xyz = gen_disordered_box(N1, N2)
            lc = init_link_cells(Lx, Nx)
            populate_link_cells(lc, xyz, Lx, Nx)
            ti = time.time()
            op = order_param_lc(types, xyz, lc, L, rc)
            tf = time.time()
            print("=== op: %.3f\nTime: %.2f s." % (op, tf - ti))
 
        if args["normal"]:
            # order
            print("Ordered box. rc: %.2f" % rc)
            types, xyz = gen_ordered_box(N1, N2)
            ti = time.time()
            op = order_param_naive(types, xyz, L, rc)
            tf = time.time()
            print("=== op: %.3f\nTime: %.2f s." % (op, tf - ti))
            if args["--save"]:
                ll.save_xyzfile("order.xyz", np.vstack((types, xyz.T)).T)
           
            # disorder
            print("Disordered box. rc: %.2f" % rc)
            types, xyz = gen_disordered_box(N1, N2)
            ti = time.time()
            op = order_param_naive(types, xyz, L, rc)
            tf = time.time()
            print("=== op: %.3f\nTime: %.2f s." % (op, tf - ti))
            if args["--save"]:
                ll.save_xyzfile("disorder.xyz", np.vstack((types, xyz.T)).T)


