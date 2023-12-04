#!/usr/bin/env python
"""
Library for manipulating boxes of atoms / particles.

Classes:
* Simulation (dt, steps, etc)
* ForceField (parameters)
* MDBox (box size, atoms, velocities, forces, bonds)

Reading from formats:
* xyz
* DL_POLY / DL_MESO
* pdb

Writing to formats:
* xyz
* DL_POLY / DL_MESO
* LAMMPS

"""
import numpy as np
import sys


class MDBox():
    def __init__(self, L=np.zeros(3)):
        if L.shape != (3,):
            sys.exit("L needs to be a vector of size 3.")
        self.L = L


    def read_config(self, fname):
        pass

    def read_xyz(self, fname):
        pass

    def read_pdb(self, fname):
        pass

    def load_xyz(self, xyz):
        pass

    def load_molecule(self, molname, bonds, bnames, Nc):
        """
        TODO: think through loading and representation of chains:
        ? How to link chains to atom numbers?
        """
        self.bonds = bonds
        self.bnames = bnames


    def gen_bin_mixt(self, N, f, types):
        """Generate binary mixture.
        * N: total number of particles
        * f: ratio of A/B particles
        * types: vector of particle types, e.g. ['A', 'B']
        """
        if len(types) != 2:
            sys.exit("There should be two types.")
        f = float(f)
        if f < 0.0 or f > 1.0:
            sys.exit("Ratio f should be between 0 and 1.")
        self.N = N
        NA = int(N * f)
        NB = N - NA
        names = [NA * types[0] + NB * types[1]]
        self.xyz = np.random.randn(N, 3) * L


    def gen_poly_solvent(self, Nc, r, Ns, types):
        """Generate polymer-solvent box
        Nc: number of chains
        r: polymerisation
        Ns: Number of solvent particles
        types: type of polymer bead, then solvent bead
        """
        self.N = Nc * r + Ns
        Np = Nc * r
        names = [Np * types[0] + Ns * types[1]]
        self.xyz = np.random.randn(self.N, 3) * L


    def gen_diblock_copolymer(self, r, NA):
        """Generate a diblock copoymer melt"""
        pass


    def gen_vel(self, kT, m=1.0):
        """Generate velocity from temperature kT"""
        self.vel = np.random.randn(self.N, 3) * kT / m


    def save_xyz(fname):
        pass




