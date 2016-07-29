#!/usr/bin/env python
"""
Collection of often used functions for
* input/output
* rotating, translating and printing data
NOW MOSTLY OBSOLETE, ONLY *_table FUNCTIONS USED

pv278@cam.ac.uk, 11/05/15
"""
import numpy as np
from numpy.matlib import repmat
from math import sin, cos


def save_xyz(coords, atom_names, filename):
    """save xyz coords into file"""
    f = open(filename,"w")
    M, N = coords.shape
    for i in range(M):
        line = str(atom_names[i])+"\t"
        for j in range(N):
            line += "%.6f" % coords[i, j] + "\t"
        line += "\n"
        f.write(line)
    f.close()
    print "Coords saved to",filename


def print_xyz(coords, atom_names):
    M, N = coords.shape
    for i in range(M):
        line = atom_names[i] + "\t"
        for j in range(N):
            line += "%.6f" % coords[i, j] + "\t"
        print line


def read_xyz(filepath):
    f = open(filepath,"r").readlines()
    A = np.array([line.split() for line in f])
    names = A[:,0]
    data = A[:,1:].astype(float)
    return names, data


def save_table(A, filepath, header=False, latex=False):
    """save table A into file, possibly with latex formatting"""
    f = open(filepath,"w")
    M,N = A.shape
    if header:
        f.write(header)
    for i in range(M):
        line = ""
        for j in range(N):
            line += str(A[i, j]) + "\t"
        if latex:
            line = "  &  ".join(line.split())
            line += "  \\\\"
        line += "\n"
        f.write(line)
    f.close()
    print "Table saved to",filepath


def print_table(A, header=""):
    if header:
        print header
    M,N = A.shape
    for i in range(M):
        line=""
        for j in range(N):
            line += str(A[i, j]) + "\t"
        print line


def read_table(filepath):
    """read summary tables
       TODO: rewrite in pandas DataFrame"""
    f = open(filepath,"r").readlines()
    A = np.array([line.rstrip("\n").split("\t") for line in f])
    return A


def get_path(Pt_dir, cluster, spin, eta=0, ext=""):
    """get full file path with eta and spin"""
    if eta != 0:
        path = Pt_dir + "/Pt_Water" + "/Pt" + cluster + "/Eta_" + str(eta) + "/S_" + str(spin) + "/Pt.out"
    else:
        path = Pt_dir + "/Pt_SP" + "/Pt" + cluster + "/S_" + str(spin) + "/Pt.out"
    if ext:
        path += "." + ext
    return path


def shift(coords, s):
    """shift coordinates by a given vector s"""
    N = len(coords)
    return coords + repmat(s,N,1)


def rotate_theta(coords, theta):
    """rotate atoms by an angle theta (in radians)"""
    N = coords.shape[0]
    Rtheta = np.array([[cos(theta),0,-sin(theta)],
                       [0,         1, 0         ],
                       [sin(theta),0, cos(theta)]])

    for i in range(N):
        coords[i,:] = np.dot(Rtheta, coords[i,:])
    return coords


def rotate_phi(coords, phi):
    """rotate atoms by angle phi (in radians)"""
    N = coords.shape[0]
    Rphi = np.array([[cos(phi),-sin(phi),0],
                     [sin(phi), cos(phi),0],
                     [0,        0,       1]])
    for i in range(N):
        coords[i,:] = np.dot(Rphi, coords[i,:])
    return coords


