#!/usr/bin/env python
"""Usage:
    kmeans_cluster.py <Nm> <frames> [--L <L> --seed <seed>]

Cluster N molecules using k-means algorithm for N/Nm centroids,
where Nm is the degree of coarse-graining (CG).
This approach does not use periodic boundary conditions.
FIX: box size.

Options:
    --L <L>        Box size [default: 29.31]
    --seed <seed>  Seed [default: 1234]

pv278@cam.ac.uk, 17/03/17
"""
import numpy as np
from scipy.cluster.vq import kmeans2, whiten
import time, glob
from docopt import docopt
from dlms_lib import read_xyzfile2, save_xyzfile


def gen_coms(xyz, clusters):
    """Generate centres of mass of the clusters"""
    Nc = len(set(clusters))
    xyz_CG = np.zeros((Nc, 3))
    cnt = 0
    for i in set(clusters):
        xyz_CG[cnt] = np.average(xyz[clusters == i], 0)
        cnt += 1
    return xyz_CG


if __name__ == "__main__":
    args = docopt(__doc__)
    seed = int(args["--seed"])
    frames = sorted(glob.glob(args["<frames>"]))
    Nf = len(frames)
    Nm = int(args["<Nm>"])
    L = float(args["--L"])

    names, xyz = read_xyzfile2(frames[0])
    N = len(xyz)
    box = np.diag(np.ones(3) * L)

    Nc = N // Nm

    print("===== K-means clustering of molecules =====")
    print("CG level: %i | Seed: %i | Frames: %i" % (Nm, seed, Nf))

    for frame in frames:
        print("\nFrame: %s" % frame)
        np.random.seed(seed)
        names, xyz = read_xyzfile2(frame)
 
        stds = np.std(xyz, 0)
        coords, pts = kmeans2(whiten(xyz), Nc)
        coords = coords * stds
        Nnz = len(set(pts))
        print("Clusters: %i | Cluster w >0 molecules: %i" % (Nc, Nnz))
        hist = np.bincount(np.bincount(pts))
        print("Number of clusters with given number of molecules:")
        print(np.c_[np.arange(len(hist)), hist].T)

        xyz_CG = gen_coms(xyz, pts)

        fname = frame.split("/")[-1].rstrip(".xyz") + "_%i.xyz" % Nm
        save_xyzfile(fname, np.c_[np.ones(len(xyz_CG)), xyz_CG])
        print("xyz file of cluster CoMs saved in %s." % fname)
 

