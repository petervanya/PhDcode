#!/usr/bin/env python
"""
Library containing a link cells class.

pv278@cam.ac.uk, 06/03/17
"""
import numpy as np
from numpy import sqrt
from numba import jit, float64, int64
import itertools
import sys, os, glob, time
from dlms_lib import save_xyzfile


@jit(float64(float64[:]), nopython=True)
def norm_numba(r):
    rn = 0.0
    for ri in r:
        rn += ri * ri
    return sqrt(rn)


class LinkCells():
    """
    Nx: number of cells in one dimension
    N: total number of cells, Nx**3
    L: box size
    dL: one cell size
    """
    def __init__(self, L=10*np.ones(3), Nx=4):
        self.Nx = Nx
        self.L = L
        self.dL = L / Nx
        self.cells = {}
        for i in range(self.Nx**3):
            self.cells[i] = []

    def __repr__(self):
        args = [(attr, getattr(self, attr)) for attr in \
                sorted( list(vars(self).keys()) )]
        return "LinkCells(%s)" % ", ".join("%s=%r" % a for a in args)


    def print(self):
        for k, v in self.cells.items():
            print(k, v)


    def gridnum2id(self, n):
        """Map 3d grid number to cell ID"""
        return int(n[0] * self.Nx**2 + n[1] * self.Nx + n[2])
    
    
    def id2gridnum(self, ID):
        """Map cell ID to 3d grid number"""
        gn = np.zeros(3).astype(int)
        gn[0] = ID // self.Nx**2
        gn[1] = (ID - gn[0] * self.Nx**2 ) // self.Nx
        gn[2] = ID - gn[0] * self.Nx**2 - gn[1] * self.Nx
        return gn


    def neighbour_cells(self, ID):
        """Find all neighbouring cells incl the given one"""
        r = self.id2gridnum(ID)
        neighs = []
        tmp = np.array([-1, 0, 1])
        for p in itertools.product(tmp, repeat=3):
            neigh = (r + p) % self.Nx
            neighs.append(neigh)
        return [self.gridnum2id(neigh) for neigh in neighs]


    def cell_pairs(self):
        """Generate list of unique pairs of neighbouring link cells"""
        cp = []
        for c in self.cells:
            cp.extend([set([c, i]) for i in self.neighbour_cells(c) if i != c])
        temp = set([frozenset(i) for i in cp])  # choose unique pairs
        self.pairs = [tuple(i) for i in temp]


    def cell_ids(self):
        """Return list of all cell IDs"""
        return self.cells.keys()


    def atom_to_cell(self, r):
        """Return cell ID for a given atom
        r: 3d array"""
        n = (r // self.dL).astype(int)
        return self.gridnum2id(n)


    def allocate_atoms(self, xyz):
        """Allocate atoms to link cells. For each atom, find out
        to which cell it belongs and append its name/number to the cell list.
        * lc: link cells"""
        N = len(xyz)
        for i in range(N):
            num = xyz[i] // self.dL % self.Nx
            self.cells[self.gridnum2id(num)].append(i)


def guess_box_size(xyz):
    """Infer box size from xyz matrix"""
    L = np.array([np.round(max(xyz[:, i]) - min(xyz[:, i]), 0) \
            for i in range(3)])
    return L


def dist_vec(xyz, box, Nx):
    """Calculate a distance vector of a xyz matrix"""

    lc = LinkCells(L=np.diag(box), Nx=Nx)
    lc.cell_pairs()
    lc.allocate_atoms(xyz)

    Nintra = sum([len(v) * (len(v) - 1) // 2 for v in lc.cells.values()])
    Ninter = sum([len(lc.cells[p[0]]) * len(lc.cells[p[1]]) for p in lc.pairs])
    dv = np.zeros(Nintra + Ninter)

    inv_box = np.linalg.pinv(box)
    G, Gn = np.zeros(3), np.zeros(3)
    dd, ddn = np.zeros(3), np.zeros(3)
    cnt = 0

    # intracell
    for cell in lc.cells.values():
        Na = len(cell)
        for i in range(Na):
            for j in range(i):
                dv[cnt] = norm_numba(xyz[cell[i]] - xyz[cell[j]])
                cnt += 1
    # intercell
    for p in lc.pairs:
       for i in lc.cells[p[0]]:
           for j in lc.cells[p[1]]:
                dd = xyz[i] - xyz[j]
                G = inv_box @ dd
                Gn = G - np.round(G)
                ddn = box @ Gn
                dv[cnt] = norm_numba(ddn)
                cnt += 1
    return dv


def test():
    """====== Test link cells ====="""
    np.random.seed(123)
    L = 10.0 * np.ones(3)
    N = 10000
    xyz = np.random.rand(N, 3) * L
    Nx = 10

    print("Testing link cells.")
    print("N = %i | L = %s | Nx = %i | Total cells = %i" % (N, L, Nx, Nx**3))
    lc = LinkCells(L, Nx)
    lc.allocate_atoms(xyz)
    
    # neighbour cells of 0
    id = 414
    print("Cell coord of id = %i (should be [4, 1, 4]):" % id, lc.id2gridnum(id))
    print("Neighbour cells of cell %i:\n" % id, lc.neighbour_cells(id))
    id = 905
    print("Cell coord of id = %i:" % id, lc.id2gridnum(id))
    print("Neighbour cells of cell %i:\n" % id, lc.neighbour_cells(id))
    Nn = len(lc.neighbour_cells(id))
    print("Total: %i (Should be 27 neighbours in total in 3D.)" % Nn)

    nc = lc.neighbour_cells(id)
    print("Neighbouring atoms:")
    na = []
    for i in nc:
        na.extend(lc.cells[i])
    print(na)

    names = np.ones((N, 1))
    mat = np.c_[names, xyz]
    save_xyzfile("cells.xyz", mat[na])


if __name__ == "__main__":
    test()
