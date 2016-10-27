#!/usr/bin/env python
"""Usage:
    density_profile.py 3d <frame> (--bead <bt> | water) [options]
    density_profile.py 2d <frame> <plane> (--bead <bt> | water) [options]
    density_profile.py test <frame> --bead <bt> [options]

Produce density profiles, bulk (3d) or slice (2d) at a given depth.
For each grid point, look for beads within 3 sigma and smear them with
f(r) ~ exp( -r^2 / 2 sigma^2)

Arguments:
    <frame>        Regex for xyz frame
    <bt>           Bead type 1..n
    <plane>        In 2d, which plane to profile on ("xy", "yz", "xz")

Options:
    --depth <d>    Where to cut the plane, relative units [default: 0.5]
    --grid <dx>    Grid spacing in DPD units [default: 0.25]
    --sigma <s>    Sigma of smearing function [default: 0.4]
    --cut <rc>     Cutoff in sigmas [default: 3]
    --nx <nx>      Number of link cells in one dim [default: 10]

pv278@cam.ac.uk, 08/09/16
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
    def __init__(self, id, number):
        self.id = id
        self.number = np.array(number, dtype=int)
        self.atoms = []

    def __repr__(self):
        args = [(attr, getattr(self, attr)) for attr in \
                sorted( list(vars(self).keys()) )]
        return "link_cell(%s)" % ", ".join("%s=%r" % a for a in args)


def create_link_cells(L, Nx):
    """
    Nx: number of cells in one dimension
    Lx: link cell size
    """
    Lx = L / Nx
    lcs = []
    cnt = 0
    for i in range(Nx):
        for j in range(Nx):
            for k in range(Nx):
                lcs.append(link_cell(cnt, [i, j, k]))
                cnt += 1
    return lcs


def allocate_atoms(lc, xyz, L, Nx):
    """Allocate atoms to link cells. For each atom, find out
    to which cell it belongs and append its name/number to the cell list.
    * lc: link cells"""
    N = len(xyz)
    Lx = L / Nx
#    print("Number of atoms to allocate: %i | LC size: %.2f" % (N, Lx))
    for i in range(N):
        num = xyz[i] // Lx % Nx
        lc[gridnum2id(num, Nx)].atoms.append(i)


def neighbour_cells(id, Nx):
    """Find all neighbouring cells incl the given one"""
    r = id2gridnum(id, Nx)
    neighs = []
    tmp = np.arange(3) - 1
    for p in itertools.product(tmp, tmp, tmp):
        neigh = (r + p) % Nx
        neighs.append(neigh)
    return [gridnum2id(neigh, Nx) for neigh in neighs]


def gridnum2id(n, Nx):
    """Map 3d grid number to cell ID"""
    return int(n[0] * Nx**2 + n[1] * Nx + n[2])


def id2gridnum(id, Nx):
    """Map cell ID to 3d grid number"""
    nx = id // (Nx**2)
    ny = (id - nx * Nx**2) // Nx
    nz = id - nx * Nx**2 - ny * Nx
    return np.array([nx, ny, nz])


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
    * r: 3d vector
    * sigma: std dev
    * rc: cutoff beyond which func = 0
    """
    nr = norm_numba(r)
    return exp(-nr**2 / (2*sigma**2)) / (2*pi*sigma**2)**(3/2) \
        if nr < rc else 0.0


def get_gridpoint(xyz, r0, lc, sigma, rc, L, Nx):
    """Generate one gridpoint. Pick local points in xyz
    and add the with the weight given by the smearing function"""
    pt = 0.0
    box = L * np.eye(3)
    inv_box = box / L**2

    id_r0 = gridnum2id((r0 // Lx), Nx)
    neigh_ids = neighbour_cells(id_r0, Nx)
    
    na = []
    for id in neigh_ids:
        na.extend(lc[id].atoms)

    for a in na:
           dr = xyz[a] - r0
           G = inv_box @ dr
           Gn = G - np.round(G)
           drn = box @ Gn
           pt += smear_func(drn, sigma, rc)
    return pt


def test_neighbour_cells():
    id = 2
    Nx = 5
    print("Testing neighbour cells.")
    print("Number of cells in one dim: %i | Cell ID = %i." % (Nx, id))
    print("Cell number (should be [0, 0, 2]):  ", id2gridnum(id, Nx))
    nc = neighbour_cells(2, Nx)
    print(nc)
    print(len(nc))


if __name__ == "__main__":
    args = docopt(__doc__)
    frame = args["<frame>"]
    if os.path.isfile(frame):
        A = ll.read_xyzfile(frame)
    else:
        print("File %s not found." % frame)
    if len(A) > 1000:   # THINK THIS THROUGH
        L = guess_box_size(A)
    else:
        L = 40.0
    sigma = float(args["--sigma"])
    rc = float(args["--cut"]) * sigma

#    Nx = int(args["--nx"])
    Nx = int(L // rc)
    Lx = L / Nx
    if rc > Lx:
        print("Warning: Smear function rc larger than link cell size.")

    dx = float(args["--grid"])
    x = np.arange(dx/2, L, dx)
    Ngrid = len(x)
    
    if args["water"]:
        xyzC = A[A[:, 0] == 3][:, 1:]
        xyzW = A[A[:, 0] == 4][:, 1:]
        lcC = create_link_cells(L, Nx)
        allocate_atoms(lcC, xyzC, L, Nx)
        lcW = create_link_cells(L, Nx)
        allocate_atoms(lcW, xyzW, L, Nx)
    if args["--bead"]:
        bead = int(args["<bt>"])
        bts = [int(b) for b in set(A[:, 0])]
        if bead not in bts:
            sys.exit("Chosen beadtype not present. Available: %s." % bts)
        xyz = A[A[:, 0] == bead][:, 1:]
        lc = create_link_cells(L, Nx)
        allocate_atoms(lc, xyz, L, Nx)


    print("===== Density profile =====")
    if args["water"]:
        print("L: %.2f | Water | Beads: %i" % (L, len(xyzC)+len(xyzW)))
    else:
        print("L: %.2f | Beadtype: %i | Beads: %i" % (L, bead, len(xyz)))
    print("Nx: %i | Total cells %i | Cell size: %.2f" % (Nx, Nx**3, Lx))
    print("Smearing sigma: %.2f | Cutoff: %.2f" % (sigma, rc))
    print("Density grid size: %.2f | Total points: %i" % (dx, Ngrid**2))


    if args["test"]:
        n = 10
        print("Atoms in link cell %i\n" % n, lc[n].atoms)
        print("Number of atoms in link cells:")
        natoms = np.array([len(lc[i].atoms) for i in range(Nx**3)])
        print(natoms)
        print("Total: %i" % sum(natoms))
        temp_xyz = xyz[lc[n].atoms]
        temp_A = np.vstack((bead*np.ones(len(temp_xyz)), temp_xyz.T)).T
        ll.save_xyzfile("test_cell.xyz", temp_A)
        print("===========")
        test_neighbour_cells()
        print("===========")


    if args["2d"]:
        planes = {"xy": (0, 1), "yz": (1, 2), "xz": (0, 2)}
        try:
            plane = planes[args["<plane>"]]
        except KeyError:
            sys.exit("Choose plane from 'xy', 'yz', 'xz'.")
        ax = set([0, 1, 2]).difference(plane)
        ax = list(ax)[0]

        d = float(args["--depth"])
        if d > 1.0 or d < 0.0:
            sys.exit("Depth should be between 0 and 1.")
        print("2d grid | Slice depth at z-coord: %.1f" % (L*d))
        rho = np.zeros((Ngrid, Ngrid))
        r0 = np.zeros(3)
        ti = time.time()
        if args["water"]:
            fname = "density_2d_water_%s_dx%.2f_x%.2f.out" % \
                    (args["<plane>"], dx, d)
            for i in range(Ngrid):
                for j in range(Ngrid):
                    r0[plane[0]] = x[i]
                    r0[plane[1]] = x[j]
                    r0[ax] = L * d
#                    r0 = np.array([x[i], x[j], L*d])
                    rho[i, j] += get_gridpoint(xyzC, r0, lcC, sigma, rc, L, Nx) / 2
                    rho[i, j] += get_gridpoint(xyzW, r0, lcW, sigma, rc, L, Nx)
        if args["--bead"]:
            fname = "density_2d_b%i_%s_dx%.2f_x%.2f.out" % \
                    (bead, args["<plane>"], dx, d)
            for i in range(Ngrid):
                for j in range(Ngrid):
                    r0 = np.array([x[i], x[j], d * L])
                    rho[i, j] = get_gridpoint(xyz, r0, lc, sigma, rc, L, Nx)
        tf = time.time()
        print("Final time: %.2f s." % (tf - ti))
        
        np.savetxt(fname, rho)
        print("2d density profile saved in %s." % fname)


    if args["3d"]:
        print("3d grid")
        rho = np.zeros((Ngrid**3, 4))
        ti = time.time()
        if args["water"]:
            fname = "density_3d_water.out"
            for i in range(Ngrid):
                for j in range(Ngrid):
                    for k in range(Ngrid):
                        r0 = np.array([x[i], x[j], x[k]])
                        gp = get_gridpoint(xyzC, r0, lcC, sigma, rc, L, Nx) / 2 + \
                             get_gridpoint(xyzW, r0, lcW, sigma, rc, L, Nx)
                        rho[i*Ngrid**2 + j*Ngrid + k] = [x[i], x[j], x[k], gp]
                print("Slice %i/%i done. %.2f s" % (i+1, Ngrid, time.time() - ti))
        if args["--bead"]:
            fname = "density_3d_b%i.out" % bead
            for i in range(Ngrid):
                for j in range(Ngrid):
                    for k in range(Ngrid):
                        r0 = np.array([x[i], x[j], x[k]])
                        rho[i*Ngrid**2 + j*Ngrid + k] = \
                            [x[i], x[j], x[k], get_gridpoint(xyz, r0, lc, sigma, rc, L, Nx)]
                print("Slice %i/%i done. %.2f s" % (i+1, Ngrid, time.time() - ti))
        tf = time.time()
        print("Final time: %.2f s." % (tf - ti))
        
        np.savetxt(fname, rho, fmt="%.6f", header="header")
        print("3d density profile saved in %s." % fname)


