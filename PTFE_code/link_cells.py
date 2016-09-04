#!/usr/bin/env python
"""
A library implementing link cells.

pv278@cam.ac.uk, 02/09/16
"""
import numpy as np
import sys
import itertools
from docopt import docopt
import lmp_lib as ll


class link_cells():
    """
    Nx: number of cells in one dimension
    N: total number of cells, Nx**3
    L: box size
    Lx: one cell size
    """
    def __init__(self, L=10.0, Nx=4):
        self.Nx = Nx
        self.N = Nx**3
        self.L = L
        self.Lx = L / Nx
        self.cells = dict()

        for i in range(self.N):   # think, if this is not too slow
            self.cells[i] = []


    def __repr__(self):
        args = [(attr, getattr(self, attr)) for attr in \
                sorted( list(vars(self).keys()) )]
        return "link_cell(%s)" % ", ".join("%s=%r" % a for a in args)


    def print_cells(self):
        """Nicely print cells and associated atoms"""
        for k, v in self.cells.items():
            print(k, v)


    def cell_ids(self):
        """Return list of all cell IDs"""
        return self.cells.keys()


    def atom_to_cell(self, r):
        """Return cell ID for a given atom
        r: 3d array"""
        n = (r // self.Lx).astype(int)
        return self.id_from_coord(n)


    def populate(self, xyz):
        """Populate cells with atom numbers"""
        for i in range(len(xyz)):
            id = self.atom_to_cell(xyz[i])
            self.cells[id].append(i)


    def cell_coord(self, id):
        """Return cell coordinates (i,j,k) from its id"""
        if id > self.N:
            sys.exit("ID must be less than the number of cells N.")
        nx = id // (self.Nx**2)
        ny = (id - nx * self.Nx**2) // self.Nx
        nz = id - nx * self.Nx**2 - ny * self.Nx
        return np.array([nx, ny, nz])


    def id_from_coord(self, n):
        """Convert cell coordinate to cell ID.
        n: 3d vector
        """
        n = np.array(n, dtype=int)
        if np.all(n > self.Nx):
            sys.exit("Each coordinate must be less than Nx.")
        return int(n[0] * self.Nx**2 + n[1] * self.Nx + n[2])


    def neighbour_cells(self, id):
        """Return IDs of neighbouring cells including
        the current cell"""
        if id not in range(self.N):
            sys.exit("ID not in the range 0..%i." % self.N)
        r = self.cell_coord(id)
        neighs = []
        tmp = np.arange(3) - 1
        for p in itertools.product(tmp, tmp, tmp):
            neigh = (r + p) % self.Nx
            neighs.append(neigh)
        return [self.id_from_coord(neigh) for neigh in neighs]


def test():
    """Test link cells on N atoms"""
    np.random.seed(123)
    L = 10.0
    N = 10000
    xyz = np.random.rand(N, 3) * L

    Nx = 10
    print("Testing link cells. ", end="")
    print("N = %i | L = %.1f | Nx = %i | cells = %i" % (N, L, Nx, Nx**3))
    lc = link_cells(L, Nx)
    # populate with atoms
    lc.populate(xyz)
    
    # neighbour cells of 0
    id = 414
    print("Cell coord of id = %i (should be [4, 1, 4]):" % id, lc.cell_coord(id))
    print("Neighbour cells of cell %i:\n" % id, lc.neighbour_cells(id))
    id = 905
    print("Cell coord of id = %i:" % id, lc.cell_coord(id))
    print("Neighbour cells of cell %i:\n" % id, lc.neighbour_cells(id))
    print("(Should be 27 neighbours in total in 3D.")

    nc = lc.neighbour_cells(id)
    print("Neighbouring atoms:")
    na = []
    for i in nc:
        na.extend(lc.cells[i])
    print(na)

    names = np.ones((N, 1))
    mat = np.hstack((names, xyz))
    ll.save_xyzfile("cells.xyz",mat[na])


if __name__ == "__main__":
    test()


