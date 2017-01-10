#!/usr/bin/env python
"""Usage: 
       dlpl_config2xyz.py <input> [--shift --shift-by <a>]
                          [--elem-names --force --vel]

Transform DL_POLY CONFIG file into xyz file.

Options:
    --shift          Shift box from centre to +ve quadrant
    --shift-by <a>   Shift by pre-defined vector
    --elem-names     Keep element names in xyz file, else use numbers
    --force          Save forces
    --vel            Save velocities

pv278@cam.ac.uk, 10/11/16
"""
import numpy as np
import sys
from docopt import docopt


def save_xyzfile(fname, names, xyz, comment="bla"):
    with open(fname, "w") as f:
        f.write(str(N) + "\n%s\n" % comment)
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

    box = np.array([line.split() for line in conf_str[1:4]]).astype(float)
    print("===== Converting DL_POLY CONFIG file to xyz =====")
    print("System size:", np.diag(box))
    conf_str = np.array(conf_str[4:])

    N = len(conf_str) // 4
    mask = np.arange(N) * 4
    names = [conf_str[i].split()[0] for i in mask]

    if not args["--elem-names"]:
        names_dict = {}
        sn = sorted(list(set(names)))
        for i in enumerate(sn):
            names_dict[i[1]] = i[0] + 1
        names = [names_dict[i] for i in names]
        print("Atom names:", names_dict)
        comment = str(names_dict)
    else:
        comment = "blah"

    xyz = np.array([[float(j) for j in conf_str[i].split()] for i in mask+1])

    if args["--shift"]:
        Ls = box[0, 0]
        print("Box size:", Ls)
        xyz += np.array([Ls, Ls, Ls])
        xyz = xyz % Ls
        print("Shifted box by Ls = %i." % Ls)
    if args["--shift-by"]:
        Ls = float(args["--shift-by"])
        xyz += np.array([Ls, Ls, Ls])
        xyz = xyz % Ls
        print("Shifted box by Ls = %i." % Ls)

    fname = infile.strip(".out") + ".xyz"
    save_xyzfile(fname, names, xyz, comment)

    if args["--vel"]:
        vel = np.array([[float(j) for j in conf_str[i].split()] \
                for i in mask+2])
        fname = infile.strip(".out") + ".vel"
        save_xyzfile(fname, names, vel)
        print("Velocity file saved in", fname)

    if args["--force"]:
        forces = np.array([[float(j) for j in conf_str[i].split()] \
                for i in mask+3])
        fname = infile.strip(".out") + ".for"
        save_xyzfile(fname, names, vel)
        print("Force file saved in", fname)


