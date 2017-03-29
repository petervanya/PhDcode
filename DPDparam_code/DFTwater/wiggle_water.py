#!/usr/bin/env python2
"""Usage:
    wiggle_water.py <infile> (--all | --bymol) [--sigma <s> --seed <seed> --save <outfile>]

[AD HOC]
Wiggle atomic coordinates by a Gaussian noise with a given sigma.

Arguments:
    <infile>         Input xyz file

Options:
    -s <outfile>     Save into xyz file
    --sigma <s>      Gaussian sigma [default: 0.1]
    --seed <seed>    Random seed [default: 123]
    --all            Wiggle all atoms randomly
    --bymol          Wiggle each water molecule separately

pv278@cam.ac.uk, 14/08/15
"""
from docopt import docopt
import numpy as np
from xyzlib import Atoms


args = docopt(__doc__)
#print args
np.random.seed(int(args["--seed"]))
sigma = float(args["--sigma"])
A = Atoms().read(args["<infile>"])

if args["--all"]:
    noise = np.random.randn(len(A), 3)*sigma
    A.coords += noise

if args["--bymol"]:
    Nmols = len(A)/3
    for i in range(Nmols):
        mol = Atoms(names=A.names[3*i:3*(i+1)], coords=A.coords[3*i:3*(i+1), :])
        tempshift = mol.coords[0]
        #mol.shift(-tempshift)   # put oxygen at [0,0,0]
        theta = np.random.randn()*sigma
        phi = np.random.randn()*sigma
        mol.rotate(theta, phi)
        #mol.shift(tempshift)
        mol.shift(np.random.randn(3, 3)*sigma)
        A.coords[3*i:3*(i+1), :] = np.copy(mol.coords)

if args["--save"]:
    A.save(args["--save"])
else:
    print A

