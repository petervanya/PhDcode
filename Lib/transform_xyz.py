#!/usr/bin/env python
"""Usage:
    transform.py <file> [-c <c>] [-a <a>] [-r <r>] [-s <s>] [-f <f>]
                      [--save <out>]

Manipulate an xyz file from command line.
Centre, shift, align or rotate an xyz structure.

Arguments:
    <file>                   Input file

Options:
    -c <c>, --centre <c>     Atom number (1-N) to centre the molecule around.
    -a <a>, --align <a>      Align s.t. given two atoms "n1 n2" lie on x-axis
    -r <r>, --rotate <r>     Rotate by angles theta and phi (in degrees)
    -s <s>, --shift <s>      Final shift of all atoms by a "x y z"
    -f <f>, --flip <f>       Flip atoms "n1 n2" by rotating whole molecule
    --save <out>             Save final coords into <out>  

pv278@cam.ac.uk, 02/10/15
"""
import numpy as np
from numpy.linalg import norm
from numpy.matlib import repmat
from scipy.linalg import expm
from math import sqrt, acos, radians
import sys
from docopt import docopt
from xyzlib import Atoms


if __name__ == "__main__":
    args = docopt(__doc__,version=1.0)
#    print args
    A = Atoms().read(args["<file>"])

    if args["--centre"]:
         n = int(args["--centre"])
         s = -A.coords[n-1]
         A.shift(s)
   
    if args["--align"]:
        n1, n2 = [int(i) for i in args["--align"].split()]
        if n1 != n2:
            A.align(n1, n2)

    if args["--flip"]:
        n1, n2 = [int(i) for i in args["--flip"].split()]
        if n1 == n2:
            sys.exit("Please choose two different atoms.")
        else:
            A.flip(n1, n2)
    
    if args["--rotate"]:
        theta, phi = [radians(float(i)) for i in args["--rotate"].split()]
        A.rotate(theta, phi)
    
    if args["--shift"]:
      s = [float(i) for i in args["--shift"].split()]
      A.shift(s)
    
    if args["--save"]:
        fname = args["--save"]
        A.save(fname)
    else:
        print(A)
   

