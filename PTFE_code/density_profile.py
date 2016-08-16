#!/usr/bin/env python
"""Usage:
    density_profile.py 3d <frame> <bt> [options]
    density_profile.py 2d <frame> <bt> [options]
    density_profile.py test <frame> <bt> [options]

#density_profile.py 2d <frame> <bt> [--depth <d> --grid <dx> --sigma <s>]

Produce density profiles, bulk (3d) or slice (2d) at a given depth.
For each grid point, look for beads within 3 sigma and add them with weight
f(r) ~ exp( -r^2 / 2 sigma^2)
TODO
====
* Add option to decide which plane to slice through
* Limit on cutoff w.r.t. link cell size

Arguments:
    <frame>             xyz frame
    <bt>                Bead type 1..n

Options:
    --depth <d>         Where to look for density cut [default: 0.5]
    --grid <dx>         Grid spacing in DPD units [default: 0.25]
    --sigma <s>         Sigma of smearing function [default: 0.4]
    --cut <rc>          Cutoff in sigmas [default: 3]
    --ncells <nlc>      Number of link cells [default: 5]

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


class link_cell(object):
    def __init__(self, id, number, rmin, rmax):
        self.id = id
        self.number = np.array(number, dtype=int)
        self.rmin = np.array(rmin)
        self.rmax = np.array(rmax)
        self.atoms = []

    def __repr__(self):
        args = [(attr, getattr(self, attr)) for attr in \
                sorted( list(vars(self).keys()) )]
        return "link_cell(%s)" % ", ".join("%s=%r" % a for a in args)


def create_link_cells(A, L, Nlcx):
    """Nlcx: number of cells in one dimension"""
    dL = L // Nlcx
    lcs = []
    cnt = 0
    for i in range(Nlcx):
        for j in range(Nlcx):
            for k in range(Nlcx):
                lcs.append(link_cell(cnt, [i, j, k], \
                                     [i*dL, j*dL, k*dL], \
                                     [(i+1)*dL, (j+1)*dL, (k+1)*dL]))
                cnt += 1
    return lcs


def allocate_atoms(A, L, linkcells, Nlcx):
    """Allocate atoms to link cells. For each atom, find out
    to which cell it belongs and append its name/number to the cell list"""
    xyz = A[:, 1:]
    N = len(xyz)
    dL = L // Nlcx
    print("Number of atoms to allocate: %i | LC size: %.2f" % (N, dL))
    for i in range(N):
        num = np.floor(xyz[i] // dL)
        linkcells[gridnum2id(num, Nlcx)].atoms.append(i)


def find_neighbour_cells(id, Nlcx):
    """Find all neighbouring cells incl the given one"""
    num = id2gridnum(id, Nlcx)
    neigh = []
    tmp = np.arange(3) - 1
    for p in itertools.product(tmp, tmp, tmp):
        n = []
        for i in range(3):
            n.append((num[i] - p[i] + Nlcx) % Nlcx)
        neigh.append(n)
    return [gridnum2id(n, Nlcx) for n in neigh]


def gridnum2id(n, Nlcx):
    """Map 3d grid number to cell ID"""
    return int(n[0] * Nlcx**2 + n[1] * Nlcx + n[2])


def id2gridnum(id, Nlcx):
    """Map cell ID to 3d grid number"""
    nx = id // (Nlcx**2)
    ny = (id - nx) // Nlcx
    nz = id - nx - ny
    return [nx, ny, nz]


def atom_in_cell(r, rmin, rmax):
    """Return True if atom with coord r is in a box
    determined by rmin and rmax"""
    truth_vec = np.zeros(3, dtype=int)
    for i in range(3):
        if r[i] > rmin[i] and r[i] < rmax[i]:
            truth_vec[i] = 1
    return np.all(truth_vec)


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
def smear_func(r, sigma, cutoff):
    """
    * r: 3d vector
    * sigma: std dev
    * cutoff: cutoff beyond which func = 0
    """
    nr = norm_numba(r)
    return exp(-nr**2 / (2*sigma**2)) / (2*pi*sigma**2)**(3/2) \
        if nr < cutoff else 0.0


def get_gridpoint2(A, r0, linkcells, sigma, cutoff, L, Nlcx):
    """Generate one gridpoint. Pick local points in A
    and add the with the weight given by the smearing function"""
    pt = 0.0
    xyz = A[:, 1:]
    box = L * np.eye(3)
    inv_box = box / L**2

    id_r0 = gridnum2id(np.floor(r0 // dL), Nlcx)
#    print("r0 =", r0, "id_r0 = ", id_r0)
#    print("lc[r0] |", linkcells[id_r0])
    neigh_ids = find_neighbour_cells(id_r0, Nlcx)
#    print("\nneigh_ids for r0:", neigh_ids)
    for id in neigh_ids:
        for a in linkcells[id].atoms:
            dr = xyz[a] - r0
            G = inv_box @ dr
            Gn = G - np.round(G)
            drn = box @ Gn
            pt += smear_func(drn, sigma, cutoff)
    return pt


def test_neighbour_cells():
    id = 2
    Nlcx = 5
    print("Testing neighbour cells.")
    print("Number of cells in one dim: %i | Cell ID = %i." % (Nlcx, id))
    print("Cell number (should be [0, 0, 2]):  ", id2gridnum(id, Nlcx))
    nc = find_neighbour_cells(2, Nlcx)
    print(nc)
    print(len(nc))


if __name__ == "__main__":
    args = docopt(__doc__)
    frame = args["<frame>"]
    if os.path.isfile(frame):
        A = ll.read_xyzfile(frame)
    else:
        print("File %s not found." % frame)
    L = guess_box_size(A)
    sigma = float(args["--sigma"])
    cutoff = float(args["--cut"]) * sigma

    Nlcx = int(args["--ncells"])
    dL = L // Nlcx
    if cutoff > dL:
        print("Warning: Smear function cutoff larger than link cell size.")

    dx = float(args["--grid"])
    x = np.arange(dx/2, L, dx)
    Ngrid = len(x)

    bead = int(args["<bt>"])
    bts = [int(b) for b in set(A[:, 0])]
    if bead not in bts:
        sys.exit("Chosen beadtype not present. Available: %s." % bts)
    A = A[A[:, 0] == bead]

    print("===== Density profile =====")
    print("L: %.2f | Beadtype: %i | Beads: %i" % (L, bead, len(A)))
    print("Total cells %i | Cell size: %.2f" % (Nlcx**3, dL))
    print("Smearing sigma: %.2f | Cutoff: %.2f" % (sigma, cutoff))
    print("Density grid size: %.2f | Total points: %i" % (dx, Ngrid**2))

    linkcells = create_link_cells(A, L, Nlcx)
    allocate_atoms(A, L, linkcells, Nlcx)

    if args["test"]:
        n = 10
        print("Atoms in link cell %i\n" % n, linkcells[n].atoms)
        print("Number of atoms in link cells:")
        natoms = np.array([len(linkcells[i].atoms) for i in range(Nlcx**3)])
        print(natoms)
        print("Total: %i" % sum(natoms))
        tempA = A[linkcells[n].atoms]
        ll.save_xyzfile("cell.xyz", tempA)
        print("===========")
        test_neighbour_cells()
        print("===========")

    if args["2d"]:
        d = float(args["--depth"])
        if d > 1.0 or d < 0.0:
            sys.exit("Depth should be between 0 and 1.")
        print("Slice depth at z-coord: %.1f" % (L*d))
        rho = np.zeros((Ngrid, Ngrid))
        tii = time.time()
        for i in range(Ngrid):
#            ti = time.time()
            for j in range(Ngrid):
                r0 = np.array([x[i], x[j], d * L])
                rho[i, j] = get_gridpoint2(A, r0, linkcells, sigma, cutoff, L, Nlcx)
            print("%i, %i / %i" % (i, j, Ngrid))
#            print(rho[i])
#            tf = time.time()
#            print("Time: %.2f s." % (tf-ti))
        tff = time.time()
        print("Final time: %.2f s." % (tff - tii))
 
        fname = "density_2d_b%i_d%.2f.out" % (bead, d)
        np.savetxt(fname, rho)
        print("2d density profile saved in %s." % fname)
     
    if args["3d"]:
        raise NotImplementedError


