#!/usr/bin/env python
"""Usage:
    density_profile2.py 3d <frame> <bt> [options]
    density_profile2.py 2d <frame> <bt> [options]
    density_profile2.py test <frame> <bt> [options]

#density_profile.py 2d <frame> <bt> [--depth <d> --grid <dx> --sigma <s>]

Produce density profiles, bulk (3d) or slice (2d) at a given depth.
For each grid point, look for beads within 3 sigma and add them with weight
f(r) ~ exp( -r^2 / 2 sigma^2)
TODO
====
* Add option to decide which plane to slice through

Arguments:
    <frame>        xyz frame
    <bt>           Bead type 1..n

Options:
    --depth <d>    Where to look for density cut [default: 0.5]
    --grid <dx>    Grid spacing in DPD units [default: 0.25]
    --sigma <s>    Sigma of smearing function [default: 0.4]
    --cut <rc>     Cutoff in sigmas [default: 3]
    --nx <nx>      Number of link cells in one dim [default: 5]

pv278@cam.ac.uk, 09/08/16
"""
import numpy as np
from numpy.linalg import norm
from numba import jit, float64, int64
from math import exp, pi, sqrt
import itertools
import sys, os
import time
from docopt import docopt
import lmp_lib as ll
from link_cells import link_cells


def guess_box_size(xyz):
    """Infer box size from xyz matrix"""
    return np.round(np.max(xyz[:, 1] - np.min(xyz[:, 1]), 0))


@jit(float64(float64[:]), nopython=True)
def norm_numba(r):
    rn = 0.0
    for ri in r:
        rn += ri*ri
    return sqrt(rn)


@jit(float64(float64[:], float64, float64), nopython=True)
def smear_func(r, sigma, rc):
    """
    * r: 3D vector
    * sigma: std dev
    * rc: rc beyond which func = 0
    """
    nr = norm_numba(r)
    return exp(-nr**2 / (2*sigma**2)) / (2*pi*sigma**2)**(3/2) \
        if nr < rc else 0.0


def get_gridpoint(xyz, r0, lc, sigma, rc):
    """Generate one gridpoint. Pick local points in xyz
    and add the with the weight given by the smearing function
    * xyz: (N, 3) matrix of bead coordinates
    * r0: 3D vector
    """
    pt = 0.0
    box = lc.L * np.eye(3)
    inv_box = box / lc.L**2
    G, Gn = np.zeros(3), np.zeros(3)
    dr, drn = np.zeros(3), np.zeros(3)
   
    id_r0 = lc.atom_to_cell(r0)
    nc = lc.neighbour_cells(id_r0)
    na = []
    for i in nc:
        na.extend(lc.cells[i])


    for a in na:
        dr = xyz[a] - r0
        G = inv_box @ dr
        Gn = G - np.round(G)
        drn = box @ Gn
        pt += smear_func(drn, sigma, rc)
#    if pt != 0.0: print("%.2f   " % pt, end="")
    return pt


if __name__ == "__main__":
    args = docopt(__doc__)
    frame = args["<frame>"]
    if os.path.isfile(frame):
        A = ll.read_xyzfile(frame)
    else:
        sys.exit("File %s not found." % frame)
    if len(A) > 1000:
        L = guess_box_size(A)
    else:
        L = 40.0
    sigma = float(args["--sigma"])
    rc = float(args["--cut"]) * sigma

    Nx = int(args["--nx"])
    if rc > (L / Nx):
        print("Warning: Smear function rc larger than link cell size.")

    dx = float(args["--grid"])
    x = np.arange(dx/2, L, dx)
    Ngrid = len(x)

    bead = int(args["<bt>"])
    bts = [int(b) for b in set(A[:, 0])]   # bead types
    if bead not in bts:
        sys.exit("Chosen beadtype not present. Available: %s." % bts)
    xyz = A[A[:, 0] == bead][:, 1:]

    print("===== Density profile =====")
    print("L: %.2f | Beadtype: %i | Beads: %i" % (L, bead, len(xyz)))
    print("Link cells. Nx: %i | Lx: %.2f" % (Nx, L / Nx))
    print("Prob func. Sigma: %.2f | rc: %.2f" % (sigma, rc))
    print("Density grid dx: %.2f | Total points: %i" % (dx, Ngrid**2))

    lc = link_cells(L, Nx)
    lc.populate(xyz)

    if args["test"]:
        n = 10
        print("Atoms in link cell %i:\n" % n, lc.cells[n])
        print("Number of atoms in link cells:")
        natoms = [len(arr) for arr in lc.cells.values()]
        print(natoms)
        print("Total: %i" % sum(natoms))

    if args["2d"]:
        d = float(args["--depth"])
        if d > 1.0 or d < 0.0:
            sys.exit("Depth should be between 0 and 1.")
        print("Slice depth at z-coord: %.1f" % (L*d))
        rho = np.zeros((Ngrid, Ngrid))
        ti = time.time()
        for i in range(Ngrid):
            for j in range(Ngrid):
                r0 = np.array([x[i], x[j], d * L])
                rho[i, j] = get_gridpoint(xyz, r0, lc, sigma, rc)
#            print("%i, %i / %i" % (i+1, j+1, Ngrid))
#            print("\nTime: %.2f s." % (tf - ti))
        tf = time.time()
        print("Final time: %.2f s." % (tf - ti))
 
        fname = "density_2d_b%i_d%.2f.out" % (bead, d)
        np.savetxt(fname, rho)
        print("2d density profile saved in %s." % fname)
     
    if args["3d"]:
        raise NotImplementedError


