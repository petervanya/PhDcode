#!/usr/bin/env python
"""
Usage: gen_O.py [-h] [-s <s> <s> <s>]
Generate oxygen coordinate

-h --help              Show this message and exit
-s s s s --shift s s s Shift coordinates by [x,y,z]

20/02/15
"""
from docopt import docopt
import numpy as np
import numpy.matlib
from math import *

def savedata(coords, filename):
    f=open(filename,"w")
    s = str("O\t")
    for j in range(3):
        s += "%.6f" % coords[j] + "\t"
    f.write(s+"\n")
    f.close()
    print "Oxygen coords saved in",filename

# ===== get input
args = docopt(__doc__)
if args["-s"]:
    coords = np.array(args["<s>"]).astype(float)
else:
    print "Default shift vector: [0,0,1] A."
    coords = np.array([0,0,1])

filename = "O.xyz"
savedata(coords, filename)
