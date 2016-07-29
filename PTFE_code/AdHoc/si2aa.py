#!/usr/bin/env python
"""Usage:
   si2aa.py <infile>

[AD HOC] Modify xyz files by converting metres to Angstroms
to conform with VMD standards.

pv278@cam.ac.uk, 11/01/16
"""
import numpy as np
from docopt import docopt
import lmp_lib as ll

args = docopt(__doc__)
A = ll.read_xyzfile(args["<infile>"])
A[:, 1:] *= 1e10
outfile = "converted.xyz"
ll.save_xyzfile(outfile, A)
print("File with units in AA saved in", outfile)
