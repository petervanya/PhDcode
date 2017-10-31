#!/usr/bin/env python
"""Usage:
    liquid_configs.py <lg> <ls>

Generate liquid configurations by taking an intersection 
of liquid-solid and liquid-gas points.

24/10/17
"""
import numpy as np
import sys
from docopt import docopt


args  = docopt(__doc__)
fname_lg = args["<lg>"]
fname_ls = args["<ls>"]
try:
    lg = open(fname_lg).readlines()
    lg = [l.strip() for l in lg]
    ls = open(fname_ls).readlines()
    ls = [l.strip() for l in ls]
except FileNotFoundError:
    sys.exit("Files not found.")

l = set(lg).intersection(set(ls))
l = list(l)
np.savetxt("configs_l.out", l, fmt="%s")
print("Intersection saved in configs_l.out.")
