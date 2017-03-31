#!/usr/bin/env python
"""Usage:
    voronoi_cluster.py <Nm> <frame> [--L <L> --seed <seed>]

Cluster N molecules into Voronoi cells of N/Nm randomly generated
points in a box, where Nm is the degree of coarse-graining 
(number of molecules in a CG particle).
Fix: box size.

Options:
    --L <L>        Box size [default: 29.31]
    --seed <seed>  Seed [default: 1234]

pv278@cam.ac.uk, 17/03/17
"""
import numpy as np
from numpy.linalg import pinv, norm
from numba import jit, float64, int64
from math import sqrt
import time, glob, sys
from docopt import docopt
from dlms_lib import read_xyzfile2, save_xyzfile


def nearest_dist(xyz_pt, pts, box):
    """Find the nearest point in pts for
    the coordinate in xyz_pt (3, 1)"""
    inv_box = pinv(box)
    Np = len(pts)
    G, Gn = np.zeros(3), np.zeros(3)
    dr, drn = np.zeros(3), np.zeros(3)

    nmin = 0
    dmin = 10.0
    for n in range(Np):
        dr = xyz_pt - pts[n]
        G = inv_box @ dr
        Gn = G - np.round(G)
        drn = box @ Gn
        d = norm(drn)
        if d < dmin:
            dmin = d
            nmin = n
    return nmin


def gen_coms(xyz, pts_d):
    """Generate centres of mass of the clusters"""
    pts_d = dict([(k, v) for k, v in pts_d.items() if len(v) != 0])
    Nc = len(pts_d)
    xyz_CG = np.zeros((Nc, 3))
    for i, v in enumerate(pts_d.values()):
        xyz_CG[i] = np.average(xyz[v], 0)
    return xyz_CG


if __name__ == "__main__":
    args = docopt(__doc__)
    seed = int(args["--seed"])
    frames = sorted(glob.glob(args["<frame>"]))
    Nf = len(frames)
    if Nf == 0: sys.exit("No frames captured.")
    Nm = int(args["<Nm>"])

    s = args["--L"].split()
    if len(s) == 1:
        L = float(s[0]) * np.ones(3)
    elif len(s) == 3:
        L = np.array(s).astype(float)
    else:
        sys.exit("<L> should have size 1 or 3.")
    box = np.diag(np.ones(3) * L)

    names, xyz = read_xyzfile2(frames[0])
    N = len(xyz)
    Np = N // Nm
    np.random.seed(seed)
    pts = np.random.rand(Np, 3) * L
    pts_d = {}

    print("==== Voronoi mesh on random points =====")
    print("CG level: %i | Seed: %i | Frames: %i" % (Nm, seed, Nf))
    print("Box: %s" % L)

    for frame in frames:
        print("\nFrame: %s" % frame)
        names, xyz = read_xyzfile2(frame)
        N = len(xyz)
        for n in range(Np):
            pts_d[n] = []
 
        ti = time.time()
        for i in range(N):
            nmin = nearest_dist(xyz[i], pts, box)
            pts_d[nmin].append(i)
        tf = time.time()
        print("Time: %.2f s." % (tf - ti))
 
        ds = [len(pts_d[i]) for i in range(Np)]
        bins = np.arange(-0.5, 10, 1.0)
        h, _ = np.histogram(ds, bins=bins)
        print(list(range(10)))
        print(h)
        eh = h @ np.arange(len(h)) / sum(h)
        ehsq = h @ np.arange(len(h))**2 / sum(h)
        print("Avg: %.2f (should be %i) | Std: %.2f" % \
                (eh, Nm, sqrt(ehsq - eh**2)))
        print("Voronoi cells with at least one molecule: %i" % sum(h[1:]))
 
        xyz_CG = gen_coms(xyz, pts_d)
 
        fname = frame.split("/")[-1].rstrip(".xyz") + "_%i.xyz" % Nm
        save_xyzfile(fname, np.c_[np.ones(len(xyz_CG)), xyz_CG])
        print("xyz file saved in %s." % fname)


