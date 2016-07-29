#!/usr/bin/env python
"""
[AD HOC] Calculate correlation of two 1D profiles.
For slabs profiles, input the ratio Lm/L capturing 
the width of the membrane Lm to the total box width L.

Usage:
    correlation.py <file1> <file2> [--slab <sl>]

Options:
    --slab <sl>         Slab with ratio [default: 1.0]

pv278@cam.ac.uk, 05/04/16
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
from docopt import docopt


def corr(a, b):
    """Correlation of vectors a, b"""
    return sum(a*b) / np.sqrt(sum(a*a) * sum(b*b))


if __name__ == "__main__":
    args = docopt(__doc__)
    a = np.loadtxt(args["<file1>"])[:, 1]
    b = np.loadtxt(args["<file2>"])[:, 1]
    sl = float(args["--slab"])
    if sl < 0.0 or sl > 1.0:
        print "Ratio <sl> must be between 0 and 1."
        sys.exit()

    if len(a) != len(b):
        print "Vectors do not have same sizes."
        sys.exit()
    
    N = len(a)
    n = int(N*sl)
    if n != N:
        ne = (N - n)/2
        a = a[ne : (N-ne)]
        b = b[ne : (N-ne)]
    else:          # get rid of weird ends
        a = a[2 : N-2]
        b = b[2 : N-2]

    # rescale so that average is 0
    a -= np.average(a)
    b -= np.average(b)

    c = corr(a, b)
    print "Correlation:", c


