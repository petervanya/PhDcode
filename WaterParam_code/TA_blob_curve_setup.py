#!/usr/bin/env python
"""Usage:
    TA_blob_curve_setup.py --pos <pos> <dmin> <dmax> <delta> [--blobnum <bn>] [--counterpoise]

[AD HOC] Generate gjf files of one water blob on triflic acid

Arguments:
    <pos>        Where on TA is the blob placed TA (number 1..6, see README.md)
    <dmin>       Minimum distance of two blobs
    <dmax>       Maximum distance of two blobs
    <delta>      Distance between successive points

Options:
    --blobnum <bn>   Which blob to use [default: 0]

pv278@cam.ac.uk, 05/10/15
"""
import numpy as np
import os, sys
from math import pi
from docopt import docopt
from xyzlib import Atoms


def save_configs(TA, blob, Drange, theta, phi, outdir, blobnum="0"):
    for d in Drange:
        xyzname = "run_blob_" + blobnum + "_d" + str(d) + ".xyz"
        xyzpath = os.path.join(outdir, xyzname)
        dist = np.array([[d, 0, 0]])
        dist = Atoms(["A"], dist)
        dist.rotate(theta, phi)
        dist = dist.coords[0].round(2)
        blob.shift(dist)
        TA_blob = triflic + blob
        blob.shift(-dist)   # STUPID TRICK
        TA_blob.save(xyzpath)


args = docopt(__doc__)
#print args
try:
    pos = args["<pos>"]
except int(pos) not in range(1, 7):
    print "Allowed position numbers are 1...6, see README."
    sys.exit()
maindir = os.path.expanduser("~/DPDcoeffs/TA_WaterBlob")
blobdir = os.path.expanduser("~/DPDcoeffs/Files/Waterblobs")
blobnum = args["--blobnum"]   # blob number used, default is 0
blobfile = os.path.join(blobdir, "waterblob_" + blobnum + ".xyz")
outdir = maindir + "/Pos_" + pos

blob = Atoms().read(blobfile)
blob.shift_com()
triflic = Atoms().read(os.path.expanduser("~/DPDcoeffs/Files/triflic.xyz"))

if args["--counterpoise"]:
    triflic.names = [i + "(fragment=1)" for i in triflic.names]
    blob.names = [i + "(fragment=2)" for i in blob.names]
    outdir = maindir + "/Pos_" + pos + "_CP"

if not os.path.exists(outdir):
    os.mkdir(outdir)

Dmin = float(args["<dmin>"])
Dmax = float(args["<dmax>"])
delta = float(args["<delta>"])
Drange = np.arange(Dmin, Dmax, delta).round(2)
print "Blob position:", pos, "\nDistances: ", Drange

if pos == "1":       # above C
    Ccoord = triflic.coords[[i for i, x in enumerate(triflic.names) if x[0] == "C"][0]]
    Drange += Ccoord[2]  # shift by z-coord
    Drange = Drange.round(2)
    theta, phi = pi/2.0, 0.0
elif pos == "2":     # below S
    theta, phi = -pi/2.0, 0.0
elif pos == "3":     # towards H
    theta, phi = 0.0, 0.0
elif pos == "4":     # opposite H
    theta, phi = 0.0, pi
elif pos == "5":     # towards O w/o H
    theta, phi = 0.0, pi/3.0
elif pos == "6":     # between O w/ H and O w/o H
    theta, phi = 0.0, 2*pi/3.0

save_configs(triflic, blob, Drange, theta, phi, outdir, blobnum)

