#!/usr/bin/env python
"""Usage:
   si2dpd.py <infile> [--rc <rc>]

[AD HOC] Modify xyz files by converting metres to DPD units
to conform with VMD standards.

Options:
    --rc <rc>      DPD units [default: 8.14e-10]
    

pv278@cam.ac.uk, 11/01/16
"""
import numpy as np
from docopt import docopt
import lmp_lib as ll

args = docopt(__doc__)
rc = float(args["--rc"])
print("rc = %.2e" % rc)

A = ll.read_xyzfile(args["<infile>"])
A[:, 1:] /= rc
outfile = "converted.xyz"
ll.save_xyzfile(outfile, A)
print("xyz frame in DPD units saved in", outfile)
