#!/usr/bin/env python
"""Usage:
    op_new.py <frames> [--rc <rc>]
    op_new.py test normal [--N <N> --f <f> --L <L> --rc <rc> --seed <s> --save]
    op_new.py test linkcells [--N <N> --f <f> --L <L> --rc <rc> --seed <s> --save]

Method 2 to compute order parameter of a binary mixture.
First, compute the total number of contacts N11, N22 and N12:
of particles of type 1, 2 within rc: OP = N12/(N11 + N12)
Then rescale w.r.t. idealised OP: where Nii = Ni*(Ni-1)/2, and Nij = Ni*Nj.

Options:
    --N <N>        Number of particles [default: 100]
    --f <f>        Fraction of A particles [default: 0.5]
    --L <L>        Box size [default: 10.0]
    --rc <rc>      Cutoff distance [default: 1.3]
    --nx <nx>      Link cells boxes in one dim [default: 5]
    --seed <s>     Random seed [default: 1234]

pv278@cam.ac.uk, 10/10/16
"""
import numpy as np
from numpy.linalg import norm
from numpy import sqrt
from numba import jit, float64, int64
import itertools
import sys, glob
import time
from docopt import docopt
from lmp_lib import read_xyzfile


def gen_ordered_box(N1, N2):
    types = np.array([1] * N1 + [2] * N2)
    f = N1 / (N1 + N2)
    xyz1 = np.random.rand(N1, 3) * [L * f, L, L]
    xyz2 = np.random.rand(N2, 3) * [L * (1 - f), L, L]
    xyz2[:, 0] += L * f
    xyz = np.vstack((xyz1, xyz2))
    return types, xyz


def gen_disordered_box(N1, N2):
    types = np.array([1] * N1 + [2] * N2)
    xyz = np.random.rand(N1+N2, 3) * L
    return types, xyz


def guess_box_size(xyz):
    """Infer box size from xyz matrix"""
    return np.round(np.max(xyz[:, 1] - np.min(xyz[:, 1]), 0))


@jit(float64(float64[:]), nopython=True)
def norm_numba(r):
    rn = 0.0
    for ri in r:
        rn += ri*ri
    return sqrt(rn)


#@jit(float64(int64[:], float64[:, :], float64, float64))
def contacts_naive(types, xyz, L, rc):
    N = len(xyz)
    box = L * np.eye(3)
    inv_box = np.linalg.pinv(box)
    dr, drn = np.zeros(3), np.zeros(3)
    G, Gn = np.zeros(3), np.zeros(3)
    N11, N22, N12 = 0, 0, 0
    
    for i in range(N):
        for j in range(i):
            dr = xyz[i] - xyz[j]
            G = inv_box @ dr
            Gn = G - np.round(G)
            drn = box @ Gn
            d = norm_numba(drn)
            if d < rc:
                if types[i] == 1 and types[j] == 1:
                    N11 += 1
                elif types[i] == 2 and types[j] == 2:
                    N22 += 1
                else:
                    N12 += 1
    return N11, N22, N12


def contacts_lc(types, xyz, lc, L, Nx, rc):
    N = len(xyz)
    box = L * np.eye(3)
    inv_box = np.linalg.pinv(box)
    dr, drn = np.zeros(3), np.zeros(3)
    G, Gn = np.zeros(3), np.zeros(3)
    Lx = L / Nx
    N11, N22, N12 = 0, 0, 0
    
    for i in range(N):
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
            d = norm_numba(drn)
            if i == j:
                continue
            if d < rc:
                if types[i] == 1 and types[j] == 1:
                    N11 += 1
                elif types[i] == 2 and types[j] == 2:
                    N22 += 1
                else:
                    N12 += 1
    return N11/2, N22/2, N12/2  # prevent double counting


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
    """Find all neighbouring cells incl the given one for a given id"""
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


def gen_op1(N11, N22, N12, rN11, rN22, rN12):
    """First method to calculate OP.
    Divide a raw OP = N12/(N11 + N22) by a reference value"""
    raw = N12 / (N11 + N22)
    ref = rN12 / (rN11 + rN22)
    print("=====\n1st method: raw = N12/(N11+N22)")
    print("raw: %.2f | ref: %.2f" % (raw, ref) )
    print("1 - raw/ref = %.3f\n=====" % (1.0 - raw / ref))


def gen_op2(N11, N22, N12, rN11, rN22, rN12):
    """Second method to calculate OP.
    Divide all contacts first by the reference contacts"""
    raw = 2 * N12 / rN12 / (N11 / rN11 + N22 / rN22)
    print("=====\n2nd method: raw = 2 * N12/rN12 / (N11/rN11 + N22/rN22)")
    print("OP: 1 - raw = %.3f\n=====" % (1 - raw))


def gen_op3(N11, N22, N12, rN11, rN22, rN12):
    """3rd method to calculate OP.
    Divide all contacts first by the reference contacts"""
    raw = N12 / (N11 * N22)
    ref = rN12 / (rN11 * rN22)
    print("=====\n3rd method: raw = N12 / (N11 * N22)")
    print("raw: %.2f | reference: %.2f" % (raw, ref))
    print("OP: 1 - raw = %.3f\n=====" % (1 - raw))


def gen_op4(N11, N22, N12, rN11, rN22, rN12):
    """4th method to calculate OP.
    Take inverse of the 2nd method."""
    raw = 2 * N12 / rN12 / (N11 / rN11 + N22 / rN22)
    print("=====\n4th method: raw = 2 * N12/rN12 / (N11/rN11 + N22/rN22)")
    print("OP: 1 / raw = %.3f\n=====" % (1.0 / raw))


if __name__ == "__main__":
    args = docopt(__doc__)
    rc = float(args["--rc"])

    if not args["test"]:
        print("===== New order parameter ====")
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
            types, xyz = A[:, 0], A[:, 1:]   # TO DO: convert types to {1, 2}
            N1, N2 = len(types == 1), len(types == 2)
            rN11 = N1*(N1-1)/2
            rN22 = N2*(N2-1)/2
            rN12 = N1*N2

            L = guess_box_size(xyz)
            Lx = L / Nx
            lc = init_link_cells(Lx, Nx)
            populate_link_cells(lc, xyz, Lx, Nx)
            N11, N22, N12 = contacts_lc(types, xyz, lc, L, Nx, rc)

#            op = 1 - 2 * N12 / rN12 / (N11 / rN11 + N22 / rN22)
            op = 1 - N12 / (N11 + N22) / (rN12 / (rN11 + rN22))
            print("%.3f   " % op, end="", flush=True)
            OP += op / Nf
        tf = time.time()
        print("\nOrder parameter: %.3f\nTime: %.2f s." % (OP, tf - ti))


    else:
        np.random.seed(int(args["--seed"]))
        L = float(args["--L"])
        f = float(args["--f"])
        N = int(args["--N"])
        N1, N2 = int(f * N), int((1-f) * N)
 
        # reference number of contacts
        rN11 = N1*(N1-1)/2
        rN22 = N2*(N2-1)/2
        rN12 = N1*N2
        
        print("===== Testing NEW order parameter script =====")
        print("Num. of AB contacts vs num. of AA plus BB contacts within rc.")
        print("N: %i | L: %.1f | f: %.2f | rc: %.2f" % (N, L, f, rc))
#        print("Ref | N11: %i | N22: %i | N12: %i" % (rN11, rN22, rN12))

        if args["linkcells"]:
            Nx = int(L // rc)
            Lx = L / Nx
 
            print("Using link cells, Nx: %i" % Nx)
            print("\nOrdered box")
            types, xyz = gen_ordered_box(N1, N2)
            lc = init_link_cells(Lx, Nx)
            populate_link_cells(lc, xyz, Lx, Nx)
            ti = time.time()
            N11, N22, N12 = contacts_lc(types, xyz, lc, L, Nx, rc)
            tf = time.time()
            gen_op1(N11, N22, N12, rN11, rN22, rN12)
            gen_op2(N11, N22, N12, rN11, rN22, rN12)
 
            print("\nDisordered box")
            types, xyz = gen_disordered_box(N1, N2)
            lc = init_link_cells(Lx, Nx)
            populate_link_cells(lc, xyz, Lx, Nx)
            ti = time.time()
            N11, N22, N12 = contacts_lc(types, xyz, lc, L, Nx, rc)
            tf = time.time()
            gen_op1(N11, N22, N12, rN11, rN22, rN12)
            gen_op2(N11, N22, N12, rN11, rN22, rN12)
 
        if args["normal"]:
            print("\nOrdered box.")
            types, xyz = gen_ordered_box(N1, N2)
            ti = time.time()
            N11, N22, N12 = contacts_naive(types, xyz, L, rc)
            tf = time.time()
            gen_op1(N11, N22, N12, rN11, rN22, rN12)
            gen_op2(N11, N22, N12, rN11, rN22, rN12)
            if args["--save"]:
                ll.save_xyzfile("order.xyz", np.vstack((types, xyz.T)).T)
 
            print("\nDisordered box.")
            types, xyz = gen_disordered_box(N1, N2)
            ti = time.time()
            N11, N22, N12 = contacts_naive(types, xyz, L, rc)
            tf = time.time()
            gen_op1(N11, N22, N12, rN11, rN22, rN12)
            gen_op2(N11, N22, N12, rN11, rN22, rN12)
            if args["--save"]:
                ll.save_xyzfile("disorder.xyz", np.vstack((types, xyz.T)).T)


