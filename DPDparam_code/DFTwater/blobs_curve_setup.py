#!/usr/bin/env python
"""Usage:
    blobs_curve_setup.py <dir> <blobnum> <dmin> <dmax> <N>

[AD HOC] Create all converged blob combinations into a gjf file
with accompanying directory structure. Main dir is ~/DPDcoeffs/Curves/

Arguments:
    <dir>        Directory in which to create to subdirectories
    <blobnum>    Maximum number to permute
    <dmin>       Minimum distance of two blobs
    <dmax>       Maximum distance of two blobs
    <N>          Number of points in between
"""
import numpy as np
from itertools import combinations
import os, sys
import pickle
from xyzlib import Atoms
from docopt import docopt

args = docopt(__doc__)
#    print args

# TESTING PICKLE LIBRARY
#    pickle.dump(args, open("args.p", "wb"))
#    args2 = pickle.load(open("args.p","rb"))
#    print args2
#
#    with open("args.txt", "w") as f:
#        f.write(str(args))
#    print eval(open("args.txt", "r").read())

maindir = os.path.abspath(args["<dir>"])
blobdir = os.path.expanduser("~/DPDcoeffs/Blobs/FinalBlobs")
Nblobs = int(args["<blobnum>"])
blobnums = range(Nblobs)
if Nblobs <= 1:
    print "Have to use more than one blob to create a two-blob system!"
    sys.exit()

Dmin = float(args["<dmin>"])
Dmax = float(args["<dmax>"])
N = int(args["<N>"])
Drange = np.linspace(Dmin, Dmax, N).round(2)
print "Range of distances: ", Drange

for c in combinations(blobnums, 2):
    combdirname = "Blobs_" + str(c[0]) + "_" + str(c[1])
    combdir = os.path.join(maindir, combdirname)
    if not os.path.isdir(combdir):  # create list of subdirs by combination
        os.makedirs(combdir)
        print "Created new directory", combdirname

    blob1 = Atoms().read(blobdir + "/waterblob_" + str(c[0]) + ".xyz")
    blob2 = Atoms().read(blobdir + "/waterblob_" + str(c[1]) + ".xyz")
    blob1.shift_com()
    blob2.shift_com()

    for d in Drange:         # in each combdir, create xyz and gjf files for each blob dist
        xyzname = "waterblobs_d" + str(d) + ".xyz"
        xyzpath = os.path.join(combdir, xyzname)
        if not os.path.exists(xyzpath):
            blob2.shift([0, 0, d])
            blob12 = blob1 + blob2
            blob2.shift([0, 0, -d])   # BAD SOLUTION, FIX THIS
            blob12.save(xyzpath)

#            header = gen_g09_header(params=g09params)  #FIX THIS
            header = "%nproc=16\n#T B3LYP/6-31G* Test\n\nSome silly text\n\n"
            output = header + str(blob12) + "\n\n"
            gjfname = "run_blobs_d" + str(d) + ".gjf"
            gjfpath = os.path.join(combdir, gjfname)
            open(gjfpath, "w").write(output)
            print "g09 input file written into", gjfpath

