#!/usr/bin/env python
"""Usage:
    config2xyz.py <input> [--shift --shift-by <a> --nafion]

Transform DL_MESO CONFIG file into xyz file.

Options:
    --shift          Shift box from centre to the +ve quadrant by guessing box size
    --shift-by <a>   Shift by pre-defined vector
    --nafion         Use "ABCWEP" bead order

pv278@cam.ac.uk, 31/05/16
"""
import numpy as np
import sys
from docopt import docopt


def save_xyzfile(fname, names, xyz):
    with open(fname, "w") as f:
        f.write(str(N) + "\nbla\n")
        for i in range(len(names)):
            f.write("%s\t%f\t%f\t%f\n" % \
                   (names[i], xyz[i, 0], xyz[i, 1], xyz[i, 2]))
    print("xyz file saved in %s." % fname)


if __name__ == "__main__":
    args = docopt(__doc__)
    infile = args["<input>"]
    try:
        conf_str = open(infile, "r").readlines()
    except FileNotFoundError:
        sys.exit("File not found: %s." % infile)

    levcfg, imcon, _ = (int(s) for s in conf_str[1].split())
    if imcon == 0:     # no box coordinates
        conf_str = np.array(conf_str[2:])
    else:              # skip over box coordinates
        box = np.array([line.split() for line in conf_str[2:5]]).astype(float)
        L = np.diag(box)
        print("System size:", L)
        conf_str = np.array(conf_str[5:])

    print("Levcfg = %i" % levcfg)
    if levcfg == 0:
        N = len(conf_str) // 2
    elif levcfg == 1:
        print("File contains velocities.")
        N = len(conf_str) // 3
    elif levcfg >= 2:
        print("File contains velocities and forces.")
        N = len(conf_str) // 4
    
    mask = np.arange(N) * (levcfg + 2)
    names = [conf_str[i].split()[0] for i in mask]
    names_dict = {}

    if args["--nafion"]:
        print("Nafion bead order 'ABCWEP' used.")
        for i, bt in enumerate("ABCWEP"):
            names_dict[bt] = i + 1
    else:
        for i in enumerate(sorted(set(names))):   # bead types to numbers
            names_dict[i[1]] = i[0] + 1

    print("Bead names: %s" % 
            ["%s: %i" % (k, v) for k, v in sorted(names_dict.items())])
    names = [names_dict[i] for i in names]

    xyz = np.array([[float(j) for j in conf_str[i].split()] for i in mask+1])

    if args["--shift"]:
        if imcon == 0:
            L = np.array([np.round((max(xyz[:, i]) - min(xyz[:, i])), 1) \
                    for i in range(3)])
        xyz += L
        xyz = xyz % L
        print("Shifted box by", L / 2)
    if args["--shift-by"]:
        s = args["--shift-by"].split()
        if len(s) == 1:
            L = float(s[0]) * np.ones(3)
        elif len(s) == 3:
            L = np.array(s).astype(float)
        else:
            sys.exit("<L> should have size 1 or 3.")
        xyz += L
        xyz = xyz % L
        print("Shifted box by", L)

    fname = infile.strip(".out") + ".xyz"
    save_xyzfile(fname, names, xyz)

    if levcfg == 1:
        vel = np.array([[float(j) for j in conf_str[i].split()] \
                for i in mask+2])
        fname = infile.strip(".out") + ".vel"
        save_xyzfile(fname, names, vel)
        print("Velocity file saved in", fname)

    if levcfg >= 2:
        vel = np.array([[float(j) for j in conf_str[i].split()] \
                for i in mask+2])
        fname = infile.strip(".out") + ".vel"
        save_xyzfile(fname, names, vel)
        print("Velocity file saved in", fname)

        forces = np.array([[float(j) for j in conf_str[i].split()] \
                for i in mask+3])
        fname = infile.strip(".out") + ".for"
        save_xyzfile(fname, names, forces)
        print("Force file saved in", fname)


