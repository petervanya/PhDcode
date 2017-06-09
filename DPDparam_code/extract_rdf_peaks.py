#!/usr/bin/env python
"""Usage: extract_rdf_peaks.py <rdffiles>

Extract position and height of the first peak of RDF
for multiple CG degrees nm.
File naming: rdf_nm3_smoothed.out

pv278@cam.ac.uk, 05/06/17
"""
import numpy as np
import glob
from docopt import docopt


args = docopt(__doc__)
rdfs = glob.glob(args["<rdffiles>"])
rdfs.sort()
if len(rdfs) == 0:
    sys.exit("No smoothed RDFs found.")

print("# nm, r_max, rdf_max")
for rdffile in rdfs:
    tmp = rdffile.split("/")[-1].rstrip(".out").split("_")
    hit = [elem for elem in tmp if "nm" in elem]
    nm = int(hit[-1].lstrip("nm"))

    A = np.loadtxt(rdffile)
    r, rdf = A[:, 0], A[:, 1]
    arg = np.argmax(rdf)
    rmax, rdfmax = r[arg], rdf[arg]
    print("%i %.3f %.3f" % (nm, rmax, rdfmax))

