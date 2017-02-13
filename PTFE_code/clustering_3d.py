#!/usr/bin/env python
"""Usage: clustering_3d.py <file> <rc> [--save]

Read water heat map, select cutoff turn into black clusters.
Use flood fill algorithm to find size and number of these clusters.

Arguments:
    <file>    2d density profile
    <rc>      Cutoff for density, usually between 0 and 3

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
    return coord[0] * N**2 + coord[1] * N + coord[2]


def id2coord(ID, N):
    i = ID // (N**2)
    j = (ID - i * N**2) // N
    k = ID - i * N**2 - j * N
    return (i, j, k)


def flood_fill(A, coord, old, new):
    """Find 1 (old), replace by 2 (new).
    * coord: (i, j, k)"""
    N = len(A)
    i, j, k = coord
    if A[coord] == new: return  # security check
    if A[coord] != old: return

    A[coord] = new
    flood_fill(A, ((i+1) % N, j, k), old, new)
    flood_fill(A, ((i-1) % N, j, k), old, new)
    flood_fill(A, (i, (j+1) % N, k), old, new)
    flood_fill(A, (i, (j-1) % N, k), old, new)
    flood_fill(A, (i, j, (k+1) % N), old, new)
    flood_fill(A, (i, j, (k-1) % N), old, new)
    return


def read_to_3d(A, rc):
    """Read matrix of coordinates into a 3d cube"""
    N = np.round(len(A)**(1/3)).astype(int)
    dx = A[1, 2] - A[0, 2]
    C = np.zeros((N, N, N)).astype(int)
    
    for line in range(len(A)):
        i = int((A[line, 0] - dx / 2) / dx)
        j = int((A[line, 1] - dx / 2) / dx)
        k = int((A[line, 2] - dx / 2) / dx)
        C[i, j, k] = 1 if A[line, 3] > rc else 0
    return C
    

if __name__ == "__main__":
    args = docopt(__doc__)
    sys.setrecursionlimit(int(1e6))
    fname = args["<file>"]
    rc = float(args["<rc>"])
    
    print("====== Clustering water in 3d =====")
    ti = time.time()
    try:
        A = np.loadtxt(fname)
    except FileNotFoundError:
        sys.exit("File not found: %s." % fname)
    tf = time.time()
    print("Reading file: %s | Time: %.2f s" % (fname, tf - ti))
    N = np.round(len(A)**(1/3)).astype(int)
    B = read_to_3d(A, rc)
    max_rho, min_rho = np.max(A[:, 3]), np.min(A[:, 3])
    dx = A[1, 2] - A[0, 2]

    Nfull = np.sum(B == 1)
    if Nfull == 0:
        sys.exit("No clusters captured, too high cutoff.")
    print("Reading file: %s" % fname)
    print("Matrix Size:", B.shape, "| Grid size dx: %.2f" % dx)
    print("Max / min density: %.2f / %.2f | cutoff: %.2f" \
            % (max_rho, min_rho, rc))
    print("Full fields: %i / %i | %.2f" % (Nfull, N**3, Nfull/N**3))

    if args["--save"]:
        figdir = "Clustering"
        if not os.path.isdir(figdir):
            os.makedirs(figdir)
        fig = plt.figure()
        ax = plt.axes()
        plt.spy(B[N//2, :, :])
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        figname = "clustering_3d_rc%.2f.png" % rc
        plt.savefig(figdir + "/" + figname, bbox_inches="tight")

    sizes = []
    nc = 0
    print("Starting clustering... ", end="", flush=True)
    ti = time.time()
    for i in range(N**3):
        start = id2coord(i, N)
        if B[start] == 1:
            flood_fill(B, start, 1, 2)
            ncurr = np.sum(B == 2)
            sizes.append((ncurr - nc) * dx**3)  # convert to DPD units
            nc = ncurr

    tf = time.time()
    print("Done. Time: %.2f s." % (tf - ti))
    sizes.sort(reverse=True)
    print(sizes)
    Pinf = np.max(sizes) / (N * dx)**3

    Nc = len(sizes)
    Savg = np.average(sizes)
    print("Number of clusters: %i | Average size: %.2f" % (Nc, Savg))
    print("Percolation cluster strength: %.3f" % Pinf)


