#!/usr/bin/env python
"""Usage: clusterng.py <file> <rc> [--save]

Read water heat map, select cutoff turn into black clusters.
Use flood fill algorithm to find size and number of these clusters.

Arguments:
    <file>    2d density profile
    <rc>      Cutoff for density between 0 and 1

Options:
    --save    Create png plot of clusters

pv278@cam.ac.uk, 27/10/16
"""
import numpy as np
import matplotlib.pyplot as plt
import time
import sys, os
from docopt import docopt


def coord2id(coord, N):
    return coord[0]*N + coord[1]


def id2coord(ID, N):
    i = ID // N
    j = ID - i * N
    return (i, j)


def flood_fill(A, coord, old, new):
    """Find 1 (old), replace by 2 (new).
    * coord: (i, j)"""
    N = len(A)
    i, j = coord
    if A[coord] == new: return  # security check
    if A[coord] != old: return

    A[coord] = new
    flood_fill(A, ((i+1) % N, j), old, new)
    flood_fill(A, ((i-1) % N, j), old, new)
    flood_fill(A, (i, (j+1) % N), old, new)
    flood_fill(A, (i, (j-1) % N), old, new)
    return


if __name__ == "__main__":
    args = docopt(__doc__)
    sys.setrecursionlimit(int(1e6))
    fname = args["<file>"]
    A = np.loadtxt(fname)
    N = len(A)
    max_rho, min_rho = np.max(A), np.min(A)
    rc = float(args["<rc>"])
    if rc < 0 or rc > 1:
        sys.exit("Enter rc between 0 and 1.")
    rc = 3.0 * rc

    B = np.zeros((N, N)).astype(int)
    B[A > rc] = 1
    Nfull = np.sum(B == 1)
    print("====== Clustering water =====")
    print("Reading file: %s" % fname)
    print("Matrix Size: (%i, %i)" % A.shape)
    print("Max / min density: %.2f / %.2f | cutoff: %.2f" \
            % (max_rho, min_rho, rc))
    print("Number of full fields: %i / %i | %.2f" % (Nfull, N**2, Nfull/N**2))

    if args["--save"]:
        figdir = "Clustering"
        if not os.path.isdir(figdir):
            os.path.makedirs(figdir)
        plt.spy(B)
        plt.axis("off")
        figname = "clustering_rc%.2f.png" % rc
        plt.savefig(figdir + "/" + figname)

    sizes = []
    nc = 0
    print("Starting clustering... ", end="", flush=True)
    ti = time.time()
    for i in range(N*N):
        start = id2coord(i, N)
        if B[start] == 1:
            flood_fill(B, start, 1, 2)
            ncurr = np.sum(B == 2)
            sizes.append(ncurr - nc)
            nc = ncurr

    tf = time.time()
    print("Done. Time: %.2f s." % (tf - ti))
    sizes.sort(reverse=True)
    print(sizes)
    Pinf = np.max(sizes) / np.sum(sizes)

    Nc = len(sizes)
    Savg = np.average(sizes)
    print("Number of clusters: %i | Average size: %.2f" % (Nc, Savg))
    print("Percolation cluster strength: %.3f" % Pinf)


