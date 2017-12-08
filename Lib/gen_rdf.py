#!/usr/bin/env python
"""Usage:
    gen_rdf.py <frames> [--bt <bt> --L <L> --bins <nb> --verbose]
    gen_rdf.py test

Compute radial distribution function from the xyz files and
for given bead types.

Arguments:
    <frames>       Regex for xyz frames

Options:
    --bt <bt>      Atom name in number [default: 1]
    --bins <nb>    Number of bins [default: 500]
    --L <L>        Box size in units of xyz files

pv278@cam.ac.uk, 08/12/17
"""
import numpy as np
from numpy import pi
import glob, sys, time
from docopt import docopt
from Fcore.f_rdf import f_rdf
from dlms_lib import read_xyzfile2


def rdf_1_type(frames, box, r, atom_types):
    """Construct RDF for beads from given xyz frames."""
    rdf = np.zeros(len(r))
    V = np.prod(np.diag(box))
    dr = r[1] - r[0]
    Nb = len(r)

    nm, xyz = read_xyzfile2(frames[0])
    N = sum(nm == atom_types[0])
    Np = N * (N - 1) // 2

    for frame in frames:
        nm, xyz = read_xyzfile2(frame)
        xyz = xyz[nm == atom_types[0]]
        rdf_raw = f_rdf2.rdf_hist(xyz, dr, box, 2 * Nb)
        rdf_raw = rdf_raw[:Nb]
        rdf += rdf_raw / (4 * pi * r**2 * dr) * V / Np
    return rdf / len(frames)


def rdf_2_types(frames, box, r, atom_types):
    rdf = np.zeros(len(r))
    V = np.prod(np.diag(box))
    dr = r[1] - r[0]
    Nb = len(r)

    nm, xyz = read_xyzfile2(frames[0])
    N1 = sum(nm == atom_types[0])
    N2 = sum(nm == atom_types[1])
    Np = N1 * N2

    for frame in frames:
        nm, xyz = read_xyzfile2(frame)
        xyz1 = xyz[nm == atom_types[0]]
        xyz2 = xyz[nm == atom_types[1]]
        rdf_raw = f_rdf2.rdf_hist_2_types(xyz1, xyz2, dr, box, 2 * Nb)
        rdf_raw = rdf_raw[:Nb]
        rdf += rdf_raw /(4 * pi * r**2 * dr) * V / Np
    return rdf / len(frames)


def guess_box_size(frames):
    """Infer the box size from the xyz frame"""
    Ntf = len(frames) // 50 + 1
    Ltrial = np.zeros((Ntf, 3))
    for i in range(Ntf):
        nm, xyz = read_xyzfile2(frames[i])
        Ltrial[i] = np.array([np.round(np.max(xyz[:, i]) - \
                np.min(xyz[:, i]), 2) for i in range(1, 4)])
    return np.max(Ltrial, 0)


def test():
    import matplotlib.pyplot as plt
    from Fcore.f_rdf2 import f_rdf2
    np.random.seed(12)
    L = 10
    N = 10000
    print("==== Testing RDF on random box ====")
    print("Beads: %i | Box size: %.1f" % (N, L))
    box = np.diag(L * np.ones(3))
    V = L**3
    xyz = np.random.rand(N, 3) * L
    nm = np.ones(N).astype(int)
    Nb = 100
    bins = np.linspace(0, L / 2, Nb + 1)
    dr = bins[1] - bins[0]
    r = dr / 2 + bins[:-1]
    Np = N * (N - 1) // 2
    print("New RDF method...")
    ti = time.time()
    rdf = f_rdf2.rdf_hist(xyz, dr, box, 2 * Nb)
    rdf = rdf[:Nb]
    print(rdf)
    rdf = rdf * V / Np / (4 * pi * r**2 * dr)
    print("Time: %.2f s." % (time.time() - ti))

    print("Old RDF method...")
    ti = time.time()
    dist_vec = f_rdf.dist_vec(xyz, box)
    rdf_old, _ = np.histogram(dist_vec, bins)
    print(rdf_old)
    rdf_old = rdf_old * V / Np / (4 * pi * r**2 * dr)
    print("Time: %.2f s." % (time.time() - ti))
    plt.plot(r, rdf, r, rdf_old, "--")
    plt.ylim(bottom=0)
    plt.show()


if __name__ == "__main__":
    args = docopt(__doc__)
    if args["<frames>"] == "test":
        test()
        sys.exit()

    frames = glob.glob(args["<frames>"])
    if len(frames) == 0:
        sys.exit("No xyz frames captured.")
    Nf = len(frames)
    N = int(open(frames[0], "r").readline())
    Nb = int(args["--bins"])

    if args["--L"]:
        s = args["--L"].split()
        if len(s) == 1:
            Ls = float(s[0]) * np.ones(3)
        elif len(s) == 3:
            Ls = np.array(s).astype(float)
        else:
            sys.exit("L should have one or three elements.")
    else:
        Ls = guess_box_size(frames)
    box = np.diag(Ls)

    bins = np.linspace(0, np.min(Ls) / 2, Nb + 1)
    dr = bins[1] - bins[0]
    r = dr / 2 + bins[:-1]
    
    atom_types = [int(i) for i in args["--bt"].split()]
    if len(atom_types) > 2:
        sys.exit("Two atom at most types allowed.")
    xyz_types = set(read_xyzfile(frames[0])[:, 0])
    if not set(atom_types).issubset(xyz_types):
        sys.exit("Requested atom types not present in the xyz files.")

    print("===== Calculating rdf =====")
    print("Atoms: %s | Bins: %i | rc %.2f | xyz frames: %i" % \
         (atom_types, Nb, rc, Nf))
    print("Box: %s" % Ls)

    ti = time.time()
    if len(atom_types) == 1:
        rdf = rdf_1_type(frames, box, r, atom_types)
        fname = "rdf_%i.out" % atom_types[0]
    elif len(atom_types) == 2:
        rdf = rdf_2_types(frames, box, r, atom_types)
        fname = "rdf_%i_%i.out" % tuple(atom_types)
    tf = time.time()
    print("Time: %.2f s." % (tf - ti))

    np.savetxt(fname, np.c_[r, rdf], fmt="%.6f")
    print("RDF saved in %s." % fname)


