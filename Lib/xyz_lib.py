#!/usr/bin/env python
"""
Collection of often used functions for
* input/output
* rotating, translating and printing data in xyz format

2nd version
pv278@cam.ac.uk, 17/06/15
"""
import numpy as np
from numpy.matlib import repmat
from numpy.linalg import norm
from scipy.linalg import expm
from math import *
import sys, os, yaml


def read_Pt_spins(infile="../Files/lowest_Pt_spins.txt"):
    """
    Read lowest spin states into a dictionary,
    either WA Goddard's or our results
    NOT USED NOW
    """
    states = [line.split(":") for line in open(infile, "r")]
    states = dict(states)
    for i in states.keys():
        states[i] = map(int, states[i].split())
    return states    


class Atoms:
    """
    Store and manipulate xyz coordinates, including:
    * atom names
    * atom coordinates
    * charge
    * spin state
    """
    def __init__(self, names=[], xyz=np.array([]), charge=0, spin=1):
        self.names = list(names)
        self.xyz = np.array(xyz, dtype=float)
        self.charge = charge
        self.spin = spin


    def __len__(self):
        """Number of atoms"""
        return len(self.xyz)


    def __repr__(self):
        """
        >>> from numpy import array
        >>> m = Atoms(["Ar"], [(0, 0, 0)])
        >>> repr(m)
        "Atoms(names=['Ar'], xyz=array([[0, 0, 0]]), charge=0, spin=1)"
        >>> str(m) == str(eval(repr(m)))
        True
        """
        kwargs = [(attr, getattr(self, attr))
                  for attr in ['names', 'xyz', 'charge', 'spin']
                  if hasattr(self, attr)]
        return "Atoms(%s)" % ', '.join('%s=%r' % kw for kw in kwargs)


    def __str__(self):
        """Print xyz onto screen"""
        if self.xyz.any():
            M, N = self.xyz.shape
            line = ""
            line += str(self.charge) + " " + str(self.spin) + "\n"
            for i in range(M):
                line += self.names[i] + "\t"
                for j in range(N):
                    line += "%.6f" % self.xyz[i, j] + "\t"
                line += "\n"
            return line.rstrip("\n")
        else:
            return ""


    def __add__(self, other):
        return Atoms(self.names + other.names, \
                     np.vstack((self.xyz, other.xyz)), \
                     self.charge, \
                     self.spin)


    def save(self, filename, vmd=False):
        """
        Save xyz xyz, charge and spin into file
        Also can save in VMD-compatible format w/ 1st line 
        with number of atoms and 2nd line a comment
        """
        with open(filename, "w") as f:
            M, N = self.xyz.shape
            if vmd:
                f.write(str(M) + "\nBlabla\n")
            else:
                f.write(str(self.charge) + " " + str(self.spin) + "\n")
            for i in range(M):
                line = str(self.names[i]) + "\t"
                for j in range(N):
                    line += "%.6f" % self.xyz[i, j] + "\t"
                line += "\n"
                f.write(line)
            f.close()
            print("Coords saved to", filename)


    def read(self, fname):
        """Read the structure from an xyz file"""
        try:
            s = open(fname, "r").readlines()
        except FileNotFoundError:
            sys.exit("File %s not found." % fname)

        if len(s) < 2:
            sys.exit("Too few lines to process.")

        if fname[-3:] != "xyz":
            sys.exit("Class", self.__class__.__name__, ": xyz file \
                    to read should have extension 'xyz'.")

        firstline = s[0].split()
        if len(firstline) == 4:     # no charge/spin line
            charge, spin = (0, 1)
            M = np.array([line.split() for line in s])
        if len(firstline) == 2:     # standard format
            charge, spin = [int(s) for s in firstline]
            M = np.array([line.split() for line in s[1:]])
        elif len(firstline) == 1:   # VMD format
            charge, spin = (0, 1)
            M = np.array([line.split() for line in s[2:]])
        names = M[:, 0]
        xyz = M[:, 1:4].astype(np.float)
        return Atoms(names, xyz, charge, spin)


    def shift(self, s):
        """Shift all atoms by a vector s"""
        assert len(s) == 3
        self.xyz = self.xyz + s    


    def rotate(self, theta=0, phi=0):
        """Rotate all atoms by angles theta and phi (radians) respectively"""
        N = len(self)
        Rtheta = np.array([[cos(theta),0,-sin(theta)],
                           [0,         1, 0         ],
                           [sin(theta),0, cos(theta)]])
        Rphi = np.array([[cos(phi),-sin(phi),0],
                         [sin(phi), cos(phi),0],
                         [0,        0,       1]])
        for i in range(N):
            self.xyz[i] = np.dot(Rtheta, self.xyz[i])
        for i in range(N):
            self.xyz[i] = np.dot(Rphi, self.xyz[i])


    def rotate_SO3(self, omega, alpha):
        """Rotate coordinates A around vector omega by angle alpha"""
        R = np.matrix([[0,        -omega[2], omega[1]],
                       [ omega[2], 0,       -omega[0]],
                       [-omega[1], omega[0], 0       ]])
        R *= alpha
        for i in range(len(self)):
            self.xyz[i] = np.dot(self.xyz[i], expm(R))


    def com(self):
        """Return centre of mass of atoms"""
#        weights = yaml.load(open(sys.path[0] + "/atomic_weights.yaml").read())
        weights_yaml = os.path.dirname(__file__) + "/atomic_weights.yaml"
        weights = yaml.load(open(weights_yaml).read())
        M = sum([weights[i] for i in self.names])
        com = np.zeros(3)
        for i in range(len(self)):
            com += weights[self.names[i]] * self.xyz[i]
        return com / M


    def shift_com(self):
        """Shift all atoms to their centre of mass"""
        self.shift(-self.com())


    def flip(self, n1, n2):
        """Flip atoms at positions n1, n2 using rotation
        n1,n2 in {1...N}"""
        N = len(self)
        if (n1 not in range(1, N+1)) or (n2 not in range(1, N+1)):
            raise ValueError
        centre = np.copy(self.xyz[n1-1])
        dist = np.copy(self.xyz[n2-1] - self.xyz[n1-1])
        s = centre + dist/2.0
        self.shift(-s)
        omega = np.array([0, dist[2], -dist[1]])
        omega = omega / norm(omega)
        self.rotate_SO3(omega, pi)
        self.shift(s)


    def align(self, n1, n2, axis=[1, 0, 0]):
        """Align the xyz so that atoms n1, n2 lie on x-axis
        n1, n2 in {1...N}"""
        N = len(self)
        if (n1 not in range(1, N+1)) or (n2 not in range(1, N+1)):
            raise ValueError
        vec12 = np.copy(self.xyz[n2-1] - self.xyz[n1-1])
        alpha = acos( np.dot(vec12, axis)/norm(vec12) )
        omega = np.cross(vec12, axis)
        omega = omega / norm(omega)
        self.rotate_SO3(omega, -alpha)


    def Rg(self):
        """Calculate the radius of gyration"""
        weights_yaml = os.path.dirname(__file__) + "/atomic_weights.yaml"
        weights = yaml.load(open(weights_yaml).read())
        com = self.com()
        rg = 0.0
        M = sum([weights[i] for i in self.names])
        for i in range(len(self)):
            rg += weights[self.names[i]] * sum((self.xyz[i] - com)**2)
        return sqrt(rg / M)



if __name__ == "__main__":
    print("===== Testing the Atoms class =====")
    xyz = np.array([[0, 0, 0], [1, 0, 2], [2, 0, 3]], dtype=float)
    atoms = Atoms(["H", "He", "Li"], xyz)
    print(atoms)
    
    s = np.arange(3)
    print("* Shifting atoms by", s)
    atoms.shift(s)
    print(atoms)

    atoms.xyz = np.copy(xyz)
    theta, phi = 90.0, 90.0
    print("* Rotating atoms by theta=%i, phi=%i:" % (theta, phi))
    atoms.rotate(radians(theta), radians(phi))
    print(atoms)

    atoms.xyz = np.copy(xyz)
    print("* Shifting atoms to its centre of mass")
    atoms.shift_com()
    print(atoms)
    
    atoms.xyz = np.copy(xyz)
    n1, n2 = 1, 2
    print("* Flipping atoms %i, %i:" % (n1, n2))
    atoms.flip(n1, n2)
    print(atoms)

    atoms.xyz = np.copy(xyz)
    n1, n2 = 1, 2
    print("* Aligning atoms %i, %i:" % (n1, n2))
    atoms.align(n1, n2)
    print(atoms)
    
    atoms.xyz = np.copy(xyz)
    print(atoms)
    print("* Radius of gyration: %.2f" % atoms.Rg())


