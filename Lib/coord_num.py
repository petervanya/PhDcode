#!/usr/bin/env python
"""Usage:
    plot_coord_num.py <rdf> <N> <L> 

Plot the coordination number from an RDF
and locate a minimum.

Arguments:
    <N>       Number of particles
    <L>       Box size, '10' or '10 8 8'

12/11/17
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
from docopt import docopt


args = docopt(__doc__)
fname = args["<rdf>"]
try:
    A = np.loadtxt(fname)
except FileNotFoundError:
    sys.exit("File %s not found." % fname)

N = int(args["<N>"])
s = args["<L>"].split()
if len(s) == 1:
    L = float(eval(s[0])) * np.ones(3)
elif len(s) == 3:
    L = np.array(s).astype(float)
else:
    sys.exit("L should have one or three elements.")

r, rdf = A[:,0], A[:,1]
rho = N / np.prod(L)
dr = r[1] - r[0]
cn = rho * 4 * np.pi * np.cumsum(r**2 * rdf * dr)

plt.plot(r, cn)
plt.xlim([0, 1.5])
plt.ylim([0, 40])
plt.grid()
figname = "cn.png"
plt.savefig(figname)
print("Coordination number plot saved in %s." % figname)

cnmin, cnmax = 1, 30
r = r[(cn > cnmin) & (cn < cnmax)]    # take only nonzero coord nums
cn = cn[(cn > cnmin) & (cn < cnmax)]
rmin = 0.5
r_temp = r[r < rmin]
cn_temp = cn[r < rmin]

if len(cn_temp) >= 2:
    cnd = np.gradient(cn_temp)
    pos = np.argmin(cnd)
    print("[0:%.1f] Derivative minimum: r = %.2f, cn = %.2f" \
            % (rmin, r_temp[pos], cn_temp[pos]))

rmin, rmax = 0.6, 1.0
r_temp = r[(r > rmin) & (r < rmax)]
cn_temp = cn[(r > rmin) & (r < rmax)]
cnd = np.gradient(cn_temp)
pos = np.argmin(cnd)
print("[%.1f:%.1f] Derivative minimum: r = %.2f, cn = %.2f" \
        % (rmin, rmax, r_temp[pos], cn_temp[pos]))


