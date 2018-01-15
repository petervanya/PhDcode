#!/usr/bin/env python
"""Usage:
    direct_corr_oz.py <rdf> <rho> [--plot]

Compute direct correlation from RDF using Fourier-Bessel transform 
for spherically symmetric functions.

31/12/17
"""
import numpy as np
from numpy import pi
from scipy.integrate import simps
from docopt import docopt


def sinc(x):
    return np.where(x > 1e-10, np.sin(x) / x, 1.0)


def forward_ft(r, hr, k):
    hk = np.zeros(k.shape)
    hk[0] = simps(hr * r**2, r) * 4*pi
    for i in range(1, len(hk)):
        hk[i] = simps(hr * sinc(k[i]*r) * r**2, r) * 4*pi
    return hk


def inverse_ft(k, hk, r):
    hr = np.zeros(r.shape)
    for i in range(len(hr)):
        hr[i] = simps(hk * sinc(k*r[i]) * k**2, k) * 4*pi / (2*pi)**3
    return hr


args = docopt(__doc__)
fname = args["<rdf>"]
try:
    A = np.loadtxt(fname)
except FileNotFoundError:
    sys.exit("File %s not found." % fname)
rho = float(args["<rho>"])
r, gr = A[:, 0], A[:, 1]
dr = r[1] - r[0]
N = len(A)
hr = gr - 1.0

k = np.array([2*pi*i / (2*N*dr) for i in range(N)])
hk = forward_ft(r, hr, k)
ck = hk / (1 + rho * hk)
cr = inverse_ft(k, ck, r)

outname = "sk.out"
np.savetxt(outname, np.c_[k, rho * hk])
print("Structure factor saved in %s." % outname)

outname = "cr_oz.out"
np.savetxt(outname, np.c_[r, cr])
print("cr saved in %s." % outname)

if args["--plot"]:
    import matplotlib.pyplot as plt
    plt.plot(r, hr, label="$h(r)$", lw=1)
    plt.plot(k/(2*pi), rho*hk, "+-", label="$\\rho h(k)$")
    plt.plot(k/(2*pi), ck, "+-", label="$c(k)$", lw=1)
    plt.plot(r, cr, label="$c(r)$")
    plt.legend()
    plt.xlabel("$r$ or $k=1/r$")
    plt.xlim([-0.1, 5])
    plt.ylim([-10, 2])
    figname = "cr_oz.png"
    plt.savefig(figname)
    print("Plot saved in %s." % figname)
    

