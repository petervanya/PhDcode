#!/usr/bin/env python
"""Usage:
    parse_config.py <input> [--shift --shift2 <a> --nafion]

Transform DL_MESO CONFIG file into xyz file.

Options:
    --shift        Shift box from centre to the +ve quadrant by inferring box size
    --shift2 <a>   Shift by pre-defined vector
    --nafion       Use "ABCWEP" bead order

pv278@cam.ac.uk, 31/05/16
"""
import numpy as np
import sys
from docopt import docopt


def save_xyzfile(fname, names, xyz):
    with open(fname, "w") as f:
        f.write(str(N) + "\nbla\n")
        for i in range(len(names)):
            f.write("%s\t%f\t%f\t%f\n" % (names[i], xyz[i, 0], xyz[i, 1], xyz[i, 2]))
    print("xyz file saved in %s" % fname)


if __name__ == "__main__":
    args = docopt(__doc__)
    infile = args["<input>"]
    try:
        conf_str = open(infile, "r").readlines()
    except FileNotFoundError:
        print("File not found:", infile)
        sys.exit()

    levcfg, imcon = (int(s) for s in conf_str[1].split())
    if imcon == 0:     # no box coordinates
        conf_str = np.array(conf_str[2:])
    else:              # skip over box coordinates
        box = np.array([line.split() for line in conf_str[2:5]]).astype(float)
        print("System size:", np.diag(box))
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
    
    mask = np.arange(N)*(levcfg + 2)
    names = [conf_str[i].split()[0] for i in mask]
    names_dict = {}

    if args["--nafion"]:
        print("Nafion bead order 'ABCWEP' used.")
        for i, bt in enumerate("ABCWEP"):
            names_dict[bt] = i+1
    else:
        for i in enumerate(set(names)):   # converting bead types to numbers
            names_dict[i[1]] = i[0]+1
    names = [names_dict[i] for i in names]

    xyz = np.array([[float(j) for j in conf_str[i].split()] for i in mask+1])

    if args["--shift"]:
        Ls = np.round((np.max(xyz) - np.min(xyz)), 1)
        xyz += np.array([Ls, Ls, Ls])
        xyz = xyz % Ls
        print("Shifted box by Ls = %i." % Ls)
    if args["--shift2"]:
        Ls = float(args["--shift2"])
        xyz += np.array([Ls, Ls, Ls])
        xyz = xyz % Ls
        print("Shifted box by Ls = %i." % Ls)

    fname = infile.strip(".out") + ".xyz"
    save_xyzfile(fname, names, xyz)

    if levcfg == 1:
        vel = np.array([[float(j) for j in conf_str[i].split()] for i in mask+2])
        fname = infile.strip(".out") + ".vel"
        save_xyzfile(fname, names, vel)
        print("Velocity file saved in", fname)

    if levcfg >= 2:
        vel = np.array([[float(j) for j in conf_str[i].split()] for i in mask+2])
        fname = infile.strip(".out") + ".vel"
        save_xyzfile(fname, names, vel)
        print("Velocity file saved in", fname)

        forces = np.array([[float(j) for j in conf_str[i].split()] for i in mask+3])
        fname = infile.strip(".out") + ".for"
        save_xyzfile(fname, names, vel)
        print("Force file saved in", fname)


