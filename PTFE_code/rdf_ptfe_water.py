#!/usr/bin/env python
"""Usage:
    rdf_ptfe_water.py <files> [--rc <rc> --L <L> --subst <s> --bins <nbins>] 

Read xyz files and compute radial distribution function
for PTFE beads (A, B, C, W, E, P). 
To get water from beads, consider:
* C: 3 molecules (SO3H + 3H2O) (bead 3)
* W: 6 molecules (bead 4)
Option to specify cutoff to consider only pairs up to this distance,
in case there are too many beads. Correct implementation of PBCs.
SO FAR ONLY WATER IMPLEMENTED, FOR OTHER BEADS CHECK rdf_ptfe.py

Arguments:
    <files>          Match the required xyz files (regex)

Options:
    --subst <s>      Substance: "water", "sulfonic", "backbone" [default: water]
    --rc <rc>        Only pairs up to this <rc> in in DPD units [default: 8]
    --bins <nbins>   Number of bins [default: 1000]
    --L <L>          Box size in DPD units [default: 40]

pv278@cam.ac.uk, 05/04/16
"""
import numpy as np
import glob
from docopt import docopt
from Fcore.f_rdf import f_rdf      # Fortran module
import lmp_lib as ll

r_DPD = 8.14e-10


def compute_rdf_water(outfile, rc, cell, bins):
    """Compute radial dist'n fcn from the xyz frame using 
    Fortran routine for pair distances and numpy binning
    """
    A = ll.read_xyzfile(outfile)
    xyz_C = A[A[:, 0] == 3][:, 1:]
    xyz_W = A[A[:, 0] == 4][:, 1:]
    NC = len(xyz_C)
    NW = len(xyz_W)
    print("NC:", NC, "| NW:", NW,"| Calculating rdf...")

    r = bins[:-1] + np.diff(bins)/2.0
    dr = bins[1] - bins[0]
    L = cell[0, 0]

    nn = int(2*NW**2 * (rc/L)**3 * 1.3)  # bulgarian const
    print("W beads, memory size:", nn*8/(1024.0)**3, "GB")
    d_W = f_rdf.dist_vec_cut(xyz_W, rc, L, cell, nn)  # rdf for W beads
    d_W = d_W[d_W != 0.0]
    print("Real size:", len(d_W)*8/(1024.0)**3, "GB")
    rdf_W, _ = np.histogram(d_W, bins)
    del d_W
    rdf_W = rdf_W /(4*np.pi*r**2 * dr) * L**3 / (NW*(NW-1)/2)

    nn = int(2*NC**2 * (rc/L)**3 * 1.3)
    print("C beads, memory size:", nn*8/(1024.0)**3, "GB")
    d_C = f_rdf.dist_vec_cut(xyz_C, rc, L, cell, nn)  # rdf for C beads
    d_C = d_C[d_C != 0.0]
    print("Real size:", len(d_C)*8/(1024.0)**3, "GB")
    rdf_C, _ = np.histogram(d_C, bins)
    del d_C
    rdf_C = rdf_C /(4*np.pi*r**2 * dr) * L**3 / (NC*(NC-1)/2)

    nn = int(2*NW**2 * (rc/L)**3 * 1.3)
    print("C and W beads, memory size:", nn*8/(1024.0)**3, "GB")
    d_CW = f_rdf.dist_vec_cut_2mat(xyz_C, xyz_W, rc, L, cell, nn)  # rdf for combined C and W beads
    d_CW = d_CW[d_CW != 0.0]
    print("Real size:", len(d_CW)*8/(1024.0)**3, "GB")
    rdf_CW, _ = np.histogram(d_CW, bins)
    del d_CW
    rdf_CW = rdf_CW /(4*np.pi*r**2 * dr) * L**3 / (NC*NW)

    norm = 6**2 + 3**2 + 6*3
    rdf = (rdf_W * 6**2 + rdf_C * 3**2 + rdf_CW * 6*3) / norm
    return rdf


def master_rdf_water(dumpfiles, rc, cell, bins):
    """Construct an rdf for water beads from given xyz frames"""
    rdf_mat = []
    for outfile in dumpfiles:
        rdf_i = compute_rdf_water(outfile, rc, cell, bins)
        rdf_mat.append(rdf_i)
        print(outfile, "done.")
    rdf_mat = np.array(rdf_mat).T
    np.savetxt("rdf_mat.out", rdf_mat)
    print("rdf matrix saved in rdf_mat.out")
    rdf = np.array(np.sum(rdf_mat, 1) / len(dumpfiles))
    return rdf


if __name__ == "__main__":
    args = docopt(__doc__)
    dumpfiles = glob.glob(args["<files>"])
    subst = args["--subst"]
    L = float(args["--L"])*r_DPD
    rc = float(args["--rc"])*r_DPD
    Nbins = int(args["--bins"])
    Nfiles = len(dumpfiles)
    N = int(open(dumpfiles[0], "r").readline())
 
    cell = L*np.eye(3)
    bins = np.linspace(0, rc, Nbins+1)
    r = bins[:-1] + np.diff(bins)/2.0
    
    if len(dumpfiles) == 0:
        raise ValueError("No xyz files captured, aborting.")
    
    print("===== Calculating rdf =====")
    print("Substance:", subst, "| Bins:", Nbins, "| Cutoff:", rc, \
          "| xyz files:", len(dumpfiles))
    
    if subst == "water":
        vals = master_rdf_water(dumpfiles, rc, cell, bins)
        fname = "rdf_water.out"
    else:
        raise NotImplementedError
#        r, vals = master_rdf(dumpfiles, rc, int(subst), Nbins)
#        fname = "rdf_" + subst + ".out"

    np.savetxt( fname, np.hstack((np.matrix(r).T, np.matrix(vals).T)) )
    print("rdf saved in", fname)



