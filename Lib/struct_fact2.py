#!/usr/bin/env python
"""Usage:
    struct_fact2.py <rdf> [--rho <rho> --sd <sd> --plot]
    struct_fact2.py <rdf> si [--rho <rho> --rc <rc> --nm <nm> --sd <sd> --plot]

Compute the structure factor from the RDF.
Sampling density (SD): k-point spacing dk = 1/(SD*L),
where L is the range of the RDF.

Options:
    --rho <rho>      Density [default: 3.0]
    --rc <rc>        DPD cutoff in nm [default: 0.814]
    --nm <nm>        CG degree [default: 6]
    --sd <sd>        Sampling density [default: 1]

pv278@cam.ac.uk, 29/12/17
"""
import numpy as np
from numpy import pi, sin
from scipy.integrate import simps
import sys, time
from docopt import docopt


def sinc(x):
    return np.where(np.abs(x) < 1e-10, 1.0, sin(x) / x)


def gen_sk(r, gr, k, rho):
    hr = gr - 1.0
    Sk = np.zeros(k.shape)
    Sk[0] = 1.0 + simps(hr * r**2, r) * 4 * pi * rho
    for i in range(1, len(Sk)):
        Sk[i] = 1.0 + simps(hr * sinc(2*pi*k[i]*r) * r**2, r) * 4 * pi * rho
    return Sk


if __name__ == "__main__":
    args = docopt(__doc__)
    try:
        A = np.loadtxt(args["<rdf>"])
    except FileNotFoundError:
        sys.exit("RDF at %s not found." % args["<rdf>"])
    rho = float(eval(args["--rho"]))
    r, gr = A[:, 0], A[:, 1]
    SD = int(args["--sd"])

    print("===== Structure factor =====")
    print("Density: %.6f" % rho)
    print("Sampling density: %i" % SD)
    if args["si"]:
        rc = float(args["--rc"])
        Nm = int(args["--nm"])
        rho = rho * Nm / rc**3
        r = r * rc
        print("Using SI units")
        print("Nm: %i | rc: %.3f nm | rho: %.2f 1/nm^3" % (Nm, rc, rho))

    N = len(r)
    dr = r[1] - r[0]
    k = np.array([i / (SD*2*N*dr) for i in range(SD*N)])

    ti = time.time()
    Sk = gen_sk(r, gr, k, rho)
    print("Time: %.2f s." % (time.time() - ti))
    k = 2 * pi * k

    outname = "sk_si.out" if args["si"] else "sk.out"
    np.savetxt(outname, np.c_[k, Sk])
    print("File saved in %s." % outname)

    if args["--plot"]:
        import matplotlib.pyplot as plt
        plt.plot(k, Sk)
        plt.ylabel("$S(k)$")
        xlbl = "$k$ (1/nm)" if args["si"] else "$k$"
        plt.xlabel(xlbl)
        plt.ylim([0.0, 3.0])
        plt.xlim([0.0, 20.0])
        plt.grid()
        figname = "sk_si.png" if args["si"] else "sk.png"
        plt.savefig(figname)
        print("S(k) plot saved in %s." % figname)
        

