#!/usr/bin/env python
"""Usage:
    smooth_rdf.py <rdffile> [--window <w> --order <o>]

Use Savitzky-Golay filter to smooth an RDF.

Options:
    --window <w>   Filter window size [default: 51]
    --order <o>    Polynomial order [default: 3]

pv278@cam.ac.uk, 07/04/17
"""
import numpy as np
from scipy.signal import savgol_filter
import sys
from docopt import docopt


if __name__ == "__main__":
    args = docopt(__doc__)
    fname = args["<rdffile>"]
    w = int(args["--window"])
    o = int(args["--order"])
    try:
        A = np.loadtxt(fname)
    except FileNotFoundError:
        sys.exit("File %s not found.")
    r, rdf = A[:, 0], A[:, 1]
    y = savgol_filter(rdf, w, o)

    outname = fname.rstrip(".out") + "_smoothed.out"
    np.savetxt(outname, np.c_[r, y])
    print("Smoothed RDF saved in %s." % outname)
    
