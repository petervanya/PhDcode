#!/usr/bin/env python
"""Usage:
    kmeans_cluster_mix.py scipy <Nm1> <Nm2> <frames> [--seed <s>]
    kmeans_cluster_mix.py sklearn <Nm1> <Nm2> <frames> [--seed <s>]

Cluster N molecules using k-means algorithm for N/Nm centroids,
where Nm is the degree of coarse-graining (number of molecules
in a CG particle).
This approach does not use periodic boundary conditions.

Options:
    --seed <s>     Seed [default: 1234]

pv278@cam.ac.uk, 21/06/17
"""
import numpy as np
from scipy.cluster.vq import kmeans2, whiten
from sklearn.cluster import KMeans
import time, glob, sys
from docopt import docopt
from dlms_lib import read_xyzfile2, save_xyzfile


if __name__ == "__main__":
    args = docopt(__doc__)
    seed = int(args["--seed"])
    np.random.seed(seed)
    frames = sorted(glob.glob(args["<frames>"]))
    Nf = len(frames)
    if Nf == 0: sys.exit("No frames captured.")
    Nm1 = int(args["<Nm1>"])
    Nm2 = int(args["<Nm2>"])

    names, xyz = read_xyzfile2(frames[0])

    print("===== K-means clustering of molecules =====")
    print("Package: %s" % "scipy" if args["scipy"] else "sklearn")
    print("CG levels: %i, %i | Seed: %i | Frames: %i" % (Nm1, Nm2, seed, Nf))

    for frame in frames:
        print("\nFrame: %s" % frame)
        names, xyz = read_xyzfile2(frame)
        xyz1, xyz2 = xyz[names == 1], xyz[names == 2]
        N1, N2 = len(xyz1), len(xyz2)
        Nc1, Nc2 = N1 // Nm1, N2 // Nm2

        ti = time.time()
        if args["scipy"]:
            stds1 = np.std(xyz1, 0)
            xyz_CG1, pts1 = kmeans2(whiten(xyz1), Nc1)
            xyz_CG1 = xyz_CG1 * stds1
            Nnz1 = len(set(pts1))
            stds2 = np.std(xyz2, 0)
            xyz_CG2, pts2 = kmeans2(whiten(xyz2), Nc2)
            xyz_CG2 = xyz_CG2 * stds2
            Nnz2 = len(set(pts2))
            print("Clusters: %i, %i | Cluster w >0 molecules: %i, %i" \
                    % (Nc1, Nc2, Nnz1, Nnz2))

        if args["sklearn"]:
            km1 = KMeans(n_clusters=Nc1).fit(xyz1)
            xyz_CG1 = km1.cluster_centers_
            pts1 = km1.labels_
            km2 = KMeans(n_clusters=Nc2).fit(xyz2)
            xyz_CG2 = km2.cluster_centers_
            pts2 = km2.labels_
        tf = time.time()
        print("Time: %.2f s." % (tf - ti))

        hist1 = np.bincount(np.bincount(pts1))
        hist2 = np.bincount(np.bincount(pts2))
        print("Number of clusters with given number of molecules:")
        print(np.c_[np.arange(len(hist1)), hist1].T)
        print(np.c_[np.arange(len(hist2)), hist2].T)

        fname = frame.split("/")[-1].rstrip(".xyz") + "_%i_%i.xyz" % (Nm1, Nm2)
        mat = np.r_[np.c_[np.ones(Nc1), xyz_CG1], \
                np.c_[2 * np.ones(Nc2), xyz_CG2]]
        save_xyzfile(fname, mat)
        print("xyz file of cluster CoMs saved in %s." % fname)
 

