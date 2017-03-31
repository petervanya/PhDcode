#!/usr/bin/env python
"""Usage:
    kmeans_cluster.py scipy <Nm> <frames> [--L <L> --seed <s>]
    kmeans_cluster.py sklearn <Nm> <frames> [--L <L> --seed <s>]

Cluster N molecules using k-means algorithm for N/Nm centroids,
where Nm is the degree of coarse-graining (number of molecules
in a CG particle).
This approach does not use periodic boundary conditions.
FIX: box size.

Options:
    --L <L>        Box size [default: 29.31]
    --seed <s>     Seed [default: 1234]

pv278@cam.ac.uk, 17/03/17
"""
import numpy as np
from scipy.cluster.vq import kmeans2, whiten
from sklearn.cluster import KMeans
import time, glob, sys
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
    Nc = N // Nm
    np.random.seed(seed)

    print("===== K-means clustering of molecules =====")
    print("Package: %s" % "scipy" if args["scipy"] else "sklearn")
    print("CG level: %i | Seed: %i | Frames: %i" % (Nm, seed, Nf))
    print("Box: %s" % L)

    for frame in frames:
        print("\nFrame: %s" % frame)
        names, xyz = read_xyzfile2(frame)

        ti = time.time()
        if args["scipy"]:
            stds = np.std(xyz, 0)
            xyz_CG, pts = kmeans2(whiten(xyz), Nc)
            xyz_CG = xyz_CG * stds
            Nnz = len(set(pts))
            print("Clusters: %i | Cluster w >0 molecules: %i" % (Nc, Nnz))

        if args["sklearn"]:
            km = KMeans(n_clusters=Nc).fit(xyz)
            xyz_CG = km.cluster_centers_
            pts = km.labels_
        tf = time.time()
        print("Time: %.2f s." % (tf - ti))

        hist = np.bincount(np.bincount(pts))
        print("Number of clusters with given number of molecules:")
        print(np.c_[np.arange(len(hist)), hist].T)

        fname = frame.split("/")[-1].rstrip(".xyz") + "_%i.xyz" % Nm
        save_xyzfile(fname, np.c_[np.ones(len(xyz_CG)), xyz_CG])
        print("xyz file of cluster CoMs saved in %s." % fname)
 

