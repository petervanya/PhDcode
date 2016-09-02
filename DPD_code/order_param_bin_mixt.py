#!/usr/bin/env python
"""Usage:
    order_param_bin_mixt.py [--N <N> --rc <rc> --save]

[AD HOC] Test order parameter of a binary mixture,
based on Goyal PhD thesis.

Options:
    --N <N>    Number of particles [default: 100]
    --rc <rc>  Cutoff distance [default: 1.3]

31/08/16
"""
import numpy as np
from numpy.linalg import norm
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


def dist_vec(types, xyz, L):
    N = len(xyz)
    Np = N * (N-1) // 2
    dv = np.zeros(Np)                   # dist vec
    tv = np.zeros((Np, 2), dtype=int)   # type vec

    cell = L * np.eye(3)
    inv_cell = np.linalg.pinv(cell)
    dr, drn = np.zeros(3), np.zeros(3)
    G, Gn = np.zeros(3), np.zeros(3)
    
    cnt = 0
    for i in range(N):
        for j in range(i):
            dr = xyz[i] - xyz[j]
            G = inv_cell @ dr
            Gn = G - np.round(G)
            drn = cell @ Gn
            dv[cnt] = norm(drn)
            tv[cnt] = [types[i], types[j]]
            cnt += 1
    return dv, tv


def order_param_naive(types, xyz, L, rc):
    N = len(xyz)
    cell = L * np.eye(3)
    inv_cell = np.linalg.pinv(cell)
    dr, drn = np.zeros(3), np.zeros(3)
    G, Gn = np.zeros(3), np.zeros(3)
    op = []
    
    for i in range(N):
        phi = 0.0
        n1, n2 = 0, 0
        for j in range(N):
            dr = xyz[i] - xyz[j]
            G = inv_cell @ dr
            Gn = G - np.round(G)
            drn = cell @ Gn
            d = norm(drn)
            if i == j:
                continue
            if d < rc:
                if types[j] == 1:
                    n1 += 1
                if types[j] == 2:
                    n2 += 2
        if n1 + n2 != 0:
            phi = (n1 - n2)**2 / (n1 + n2)**2
            op.append(phi)
#        else:
#            print("WARNING: no neighbours within given cutoff.")
#    print(len(op))
    return sum(op) / len(op)


if __name__ == "__main__":
    args = docopt(__doc__)
    np.random.seed(1234)
    L = 10.0
    f = 0.5
    N = int(args["--N"])
    N1, N2 = int(f * N), int((1-f) * N)
    rc = float(args["--rc"])

    # order
    print("N: %i | L: %.1f | f: %.2f" % (N, L, f))
    # order
    print("Testing ordered box. rc: %.2f" % rc)
    ti = time.time()
    types, xyz = gen_ordered_box(N1, N2)
    op = order_param_naive(types, xyz, L, rc)
    print("=== op: %.2f" % op)
    tf = time.time()
    print("Time: %.2f s." % (tf - ti))
    
    if args["--save"]:
        ll.save_xyzfile("order.xyz", np.vstack((types, xyz.T)).T)
   
    # disorder
    print("Testing disordered box. rc: %.2f" % rc)
    ti = time.time()
    types, xyz = gen_disordered_box(N1, N2)
    op = order_param_naive(types, xyz, L, rc)
    print("=== op: %.2f" % op)
    tf = time.time()
    print("Time: %.2f s." % (tf - ti))
   
    if args["--save"]:
        ll.save_xyzfile("disorder.xyz", np.vstack((types, xyz.T)).T)


