#!/usr/bin/env python
"""Usage:
    direct_corr.py gmff <rdf> <i> <a> [--beta <b> --cut <c> --plot]
    direct_corr.py lj <rdf> <eps> <sigma> [--beta <b> --cut <c> --plot]

Produce the direct correlation function c(r) from the RDF
and potential using HNC or PY.
<i> is the potential exponent and <a> the interaction parameter.

Options:
    --beta <b>     Inverse temperature [default: 1.0]
    --cut <c>      Cutoff for low rdf values to avoid inf [default: 0.01]

10/12/17
"""
import numpy as np
from docopt import docopt


def pot_gmff(r, a, i):
    return a * (1 - r)**(i+1) / (i+1) if r <= 1 else 0.0


def pot_lj(r, eps, s):
    return 4 * eps * ((s / r)**12 - (s / r)**6)


args = docopt(__doc__)
fname = args["<rdf>"]
A = np.loadtxt(fname)
r, gr = A[:, 0], A[:, 1]

beta = float(args["--beta"])
if args["gmff"]:
    a = float(args["<a>"])
    ii = int(args["<i>"])
elif args["lj"]:
    eps = float(args["<eps>"])
    sigma = float(args["<sigma>"])

cut = float(args["--cut"])
r = r[gr > cut]
gr = gr[gr > cut]

hr = gr - 1.0
if args["gmff"]:
    vr = np.array([pot_gmff(ri, a, ii) for ri in r])
elif args["lj"]:
    vr = np.array([pot_lj(ri, eps, sigma) for ri in r])

cr = hr - beta * vr - np.log(gr)  # HNC
#cr = (1.0 - np.exp(beta * vr)) * gr # PY

if args["--plot"]:
    import matplotlib.pyplot as plt
    plt.plot(r, cr, r, hr)
    plt.show()

tmp = fname.split("/")
outname = "cr_hnc.out"
#outname = "/".join(tmp[:-1]) + "/cr.out"
np.savetxt(outname, np.c_[r, cr])
print("File saved in %s." % outname)


