#!/usr/bin/env python
"""Usage:
    xyz2config.py <xyzfile> [--box <L>]

Convert an xyz file to DL_MESO config file.
Only positions are considered.

Options:
    --box <L>      Box size [default: 10]

pv278@cam.ac.uk, 08/09/17
"""
import numpy as np
import string
from collections import OrderedDict
from docopt import docopt
from dlms_lib import read_xyzfile2, save_config


if __name__ == "__main__":
    args = docopt(__doc__)
    fname = args["<xyzfile>"]
    s = args["--box"].split()
    if len(s) == 1:
        box = float(s[0]) * np.ones(3)
    elif len(s) == 3:
        box = np.array(s).astype(float)
    else:
        sys.exit("<L> should be a vector size 1 or 3.")

    nm, xyz = read_xyzfile2(fname)
    names = set(nm)
    letters = string.ascii_uppercase
    nm2let = OrderedDict()
    for i in names:
        nm2let[i] = letters[i-1]
    print("Conversion to letters: %s" % nm2let)
    nml = [nm2let[i] for i in nm]

    outname = "CONFIG"
    save_config(outname, nml, xyz, box, imcon=0)
