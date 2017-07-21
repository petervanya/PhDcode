#!/usr/bin/env python
"""Usage:
    predict_as_mo.py <paramfile> <rdf> [--algo <algo> --cut <rc>]

Generate interaction params using maximum-overlap method
with an RDF of at given CG degree.
Optimise also w.r.t. horizontal position of the first peak,
which should be rescaled by c.a. 0.1.

Options:
    --algo <algo>    Minimisation algorithm [default: Nelder-Mead]
    --cut <rc>       Fit rdf only above rc (in AA) [default: 0.0]

28/06/17
"""
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import sys, time
from docopt import docopt


def parse_nm(s):
    tmp = s.split("/")[-1].rstrip(".out").split("_")
    hit = [elem for elem in tmp if "nm" in elem]
    return int(hit[-1].lstrip("nm"))


def par_cubic(x, a, b, c, d):
    return a * x**3 + b * x**2 + c * x + d

def par_quart(x, a, b, c, d, e):
    return a * x**4 + b * x**3 + c * x**2 + d * x + e


def func_rdf(x, params, a):
    c = np.zeros(8)
    for i in range(8):
        if params.shape[1] == 4:
            c[i] = par_cubic(a, *params[i-1])
        if params.shape[1] == 5:
            c[i] = par_quart(a, *params[i-1])

    rect = (np.tanh(c[6] * (x - c[7])) + 1.0) / 2
    f = (c[1] * np.exp(- c[2] * (x - c[3])**2) * \
            np.sin(c[4] * (x - c[5])) + 1.0) * rect
    return f


def err(a):
    """
    a[0]: interaction parameters
    a[1]: horizontal scale
    """
    f = np.array([func_rdf(a[1] * ri, params, a[0]) for ri in r])
    return np.sum((rdf - f)**2)


args = docopt(__doc__)
try:
    params = np.loadtxt(args["<paramfile>"])
except FileNotFoundError:
    sys.exit("File %s not found." % args["<paramfile>"])
rdffile = args["<rdf>"]
algo = args["--algo"]

Nm = parse_nm(rdffile)
A = np.loadtxt(rdffile)
r, rdf = A[:, 0], A[:, 1]
rc = float(args["--cut"])
rdf = rdf[r > rc]
r = r[r > rc]
print("===== Generating a using maximum overlap ====")
print("RDF file: %s\nCG degree: %i" % (rdffile, Nm))
print("Minimisation algo: %s | Bottom cutoff: %.2f AA" % (algo, rc))

a0 = [50.0, 0.3]
bnds = ((0.0, 200.0), (0.0, 0.5))
ti = time.time()

res = minimize(err, a0, method=algo, bounds=bnds)

tf = time.time()
a, scale = res.x
print(res)
print("Time: %.3f s." % (tf - ti))
print("Nm: %i | a: %.3f | horizontal scale: %.3f" % (Nm, a, scale))

ffit = np.array([func_rdf(scale * ri, params, a) for ri in r])
plt.plot(r, rdf, "--", r, ffit, "-")
plt.xlim([0, 20])
plt.ylim([0, 2])
plt.grid()
figname = "rdf_mo_nm%i.png" % Nm
plt.savefig(figname)
print("Figure saved in %s." % figname)


