#!/usr/bin/env python
"""Usage:
    test_f_rdf.py <N> [--L <L> --plot]

[AD HOC] Test the f_rdf.f90 code

Options:
    --L <L>             Box size [default: 1.0]
07/04/16
"""
import numpy as np
import matplotlib.pyplot as plt
import time
from docopt import docopt
from f_rdf import f_rdf


args = docopt(__doc__)
N = int(args["<N>"])
L = float(args["--L"])
rc = 0.3
A = np.random.rand(N, 3)
nn = int(2*N**2 * (rc/L)**3 * 1.1)

print "Particles:", N, "| Cutoff:", rc, "| Box size:", L
print "Pairs:", N*(N-1)/2, "| Predicted length of dist vec:", nn
cell = L*np.eye(3)

t1 = time.time()
f = f_rdf.pair_dist_arr_cut(A, rc, L, cell, nn)#, N*(N-1)/2)
t2 = time.time()
f = f[f != 0.0]
print "Real length of dist vec:", len(f)
print "Time:", t2-t1
h, r = np.histogram(f, bins=50)
r = r[:-1] + np.diff(r)/2.0
dr = r[1] - r[0]
print h

if args["--plot"]:
    h = h * L**3/(N*(N-1)/2) / (4 * np.pi * r**2 * dr)
    plt.plot(r, h)
    plt.ylim([0.0, 2.0])
    plt.show()
