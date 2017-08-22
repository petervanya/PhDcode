#!/usr/bin/env python
"""
Usage:
    struct_fact.py <rdf> [--rho <rho> --L <L> --plot --f <f>]
    struct_fact.py <rdf> si [--rho <rho> --L <L> --rc <rc> --nm <nm> --plot --f <f>]

Compute the structure factor from the RDF.

Options:
    --L <L>        Box size [default: 10]
    --rho <rho>    DPD density [default: 3.0]
    --rc <rc>      DPD cutoff [default: 0.814]
    --nm <nm>      CG degree [default: 6]
    --f <f>        Sampling freq, multiple of N-points to use for k's [default: 1]

pv278@cam.ac.uk, 21/08/17
"""
import numpy as np
from numpy import pi, sin
import sys
from docopt import docopt


def gen_sk(r, gr, k, rho=3.0):
    dr = r[1] - r[0]
    print(rho)
    print(dr)
    N = len(r)
    Sk = [1 + 4 * pi * rho * sum([r[j]**2 * (gr[j] - 1.0) * \
            sin(ki * r[j]) / (ki * r[j]) * dr for j in range(N)]) for ki in k]
    return np.array(Sk)


if __name__ == "__main__":
    args = docopt(__doc__)
    try:
        A = np.loadtxt(args["<rdf>"])
    except FileNotFoundError:
        sys.exit("RDF at %s not found." % args["<rdf>"])
    rho = float(args["--rho"])
    L = float(args["--L"])
    f = int(args["--f"])
    r, gr = A[:, 0], A[:, 1]

    print("===== Structure factor =====")
    print("Density: %.2f" % rho)
    if args["si"]:
        rc = float(args["--rc"])
        Nm = int(args["--nm"])
        rho = rho * Nm / rc**3
        r = r * rc
        print("Using SI units")
        print("Nm: %i | rc: %.2f nm | rho: %.2f 1/nm^3" % (Nm, rc, rho))

    N = len(r)
    dr = r[1] - r[0]
    kmin = 2 * pi / max(r)
    kmax = 2 * pi / dr
#    k = np.linspace(0, kmax, f*N+1)[1:]
    k = np.linspace(0, max(r), f*N+1)[1:]
    Sk = gen_sk(r, gr, k, rho)
    outname = "sk_si.out" if args["si"] else "sk.out"
    np.savetxt(outname, np.c_[k, Sk])
    print("File saved in %s." % outname)

    if args["--plot"]:
        import matplotlib.pyplot as plt
        plt.plot(k, Sk)
        plt.ylabel("$S(k)$")
        xlbl = "$k$ (1/nm)" if args["si"] else "$k$"
        plt.xlabel(xlbl)
        plt.ylim([0.0, 1.5])
        plt.xlim([0.0, 20.0])
        plt.grid()
        figname = "sk_si.png" if args["si"] else "sk.png"
        plt.savefig(figname)
        print("S(k) plot saved in %s." % figname)
        

