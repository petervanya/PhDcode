#!/usr/bin/env python
"""Usage:
    castep_mnpl.py create --pos <posfile> --lat <latfile> [--kpoints <kpts>] [--save <fname>]

Collection of functions to generate a Castep input file
TO DO:
* add symmetry
* express xyz file coordinates in the lattice units
* think of other possible functions manipulate Castep files (e.g. extracting output)

Arguments:
    --pos <posfile>    Atomic positions file
    --lat <latfile>    Lattice vectors file

Options:
    --kpoints <kpts>   K points [default: 8 8 8]
    --save <fname>     Save into file 

pv278@cam.ac.uk, 12/11/15
"""
import numpy as np
import sys, os, glob
import yaml
from docopt import docopt
from xyzlib import Atoms   # using class Atoms to manipulate xyz matrices

def print_lat(lat):
    """Produce string from lattice vectors"""
    s = ""
    # FILL
    return s


def print_pos(pos):
    """Produce string from xyz matrix in the form of Atoms class"""
    # TO DO: EXCLUDE FIRST LINE
    s = "%BLOCK POSITIONS_FRAC\n" + str(pos) + "\n%ENDBLOCK POSITIONS_FRAC"
    return s


def print_masses(names, weights):
    """Generate a list of atom/mass pairs from the simulated atoms"""
    unique = set(names)   # all simulated atom types, e.g. C, H
    s = "%BLOCK SPECIES_MASS\n"
    s += "\n".join(["%s\t%f" % (atom, weights[atom]) for atom in unique])
    s += "\n%ENDBLOCK SPECIES_MASS"
    return s


if __name__ == "__main__":
    args = docopt(__doc__)
#    print args
    if args["create"]:
        # TO DO: USE try...except STRUCT TO OPEN FILES
        kpts = np.array(args["--kpoints"].split(), dtype=int)
        lat = np.loadtxt(args["<latfile>"])
        pos = Atoms().read(args["<posfile>"])
 
        lat_txt = print_lat(lat)
        pos_txt = print_pos(pos)
        kpts_txt = "KPOINTS_MP_GRID %i %i %i" % (kpts[0], kpts[1], kpts[2])
        weights = yaml.load(open("atomic_numbers.txt"))
        mass_txt = print_masses(pos.names, weights)
 
        final_string = lat_txt + "\n\n" + \
                       pos_txt + "\n\n" + \
                       kpts_txt + "\n\n" + \
                       mass_txt + "\n\n"
 
        if args["--save"]:
            fname = args["--save"]
            open(fname, "w").write(final_string)
        else:
            print final_string


