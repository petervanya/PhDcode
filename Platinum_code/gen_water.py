#!/usr/bin/env python
"""
Usage:
    gen_water [-t <t>] [-f <phi>] [-s <s>] [-p <p>] [--save <file>]

Generate xyz water coordinates in xz plane with the possibility of rotation.

Options:
    -h, --help              Show this message and exit
    -s <s>,--shift <s>      Shift of O above Pt atom in "x y z"
    -p <p>,--posPt <p>      Position of O atom in Pt lattice vectors in "x y z"
    -f <phi>,--phi=<phi>    Rotation around z-axis (phi) in degrees (done before theta)
    -t <t>,--theta=<t>      Rotation around y-axis (theta) in degrees
    --save <file>           Save into specified file

pv278@cam.ac.uk, 20/02/15
"""
import numpy as np
from docopt import docopt
from math import sin, cos, sqrt
from iolib import save_xyz, print_xyz

def init_water():
    """Initialise water molecule
       Already optimised using 6-311G**"""
    coords = np.zeros((3,3))
    coords[1,0] += 0.757009                   # Before: l_OH*sin(alpha/2)
    coords[2,0] += -0.757009                  # Before: -l_OH*sin(alpha/2)
    coords[1:3,2] += 0.593565                 # Before: l_OH*cos(alpha/2)
    return coords

def shift(coords, s):
    """shift H2O atoms"""
    print "Shift by",s
    return coords + s

def shift_Pt(coords, vectPt):
    """shift H2O atoms in Pt lattice vectors"""
    aPt = 2.775
    vPt = aPt*sqrt(3.0)/2
    hPt = aPt*sqrt(2.0/3)
    basis = np.array([[aPt, aPt/2,  aPt/2],
                      [0.0, vPt,    vPt/3],
                      [0.0, 0.0,    hPt  ]])
    s = np.dot(basis, vectPt)
    print "Position on Pt lattice:", s
    return coords + s

def rotate_theta(coords, theta):
    Rtheta = np.array([[cos(theta),0,-sin(theta)],
                       [0,         1, 0         ],
                       [sin(theta),0, cos(theta)]])

    for i in range(3):
        coords[i,:] = np.dot(Rtheta, coords[i,:])
    return coords

def rotate_phi(coords, phi):
    Rphi = np.array([[cos(phi),-sin(phi),0],
                     [sin(phi), cos(phi),0],
                     [0,        0,       1]])
    for i in range(3):
        coords[i,:] = np.dot(Rphi, coords[i,:])
    return coords


if __name__ == "__main__":
    args = docopt(__doc__, version=1.0)
#    print args
    
    coords = init_water()

    if args["--phi"]:
        phi = radians(float(args["--phi"]))
        coords = rotate_phi(coords, phi)
        print "Rotated by phi =", phi

    if args["--theta"]:
        theta = radians(float(args["--theta"]))
        coords = rotate_theta(coords, theta)
        print "Rotated by theta =", theta

    if args["--posPt"]:
        vecPt = np.array(args["--posPt"].split()).astype(float)
        coords = shift_Pt(coords, vecPt)

    if args["--shift"]:
        s = np.array(args["--shift"].split()).astype(float)
        coords = shift(coords, s)
    
    names = ["O","H","H"]
    if args["--save"]:
        filename = args["--save"]
        save_xyz(coords, names, filename)
    else:
        print_xyz(coords, names)
