#!/usr/bin/env python
"""Usage:
    parse_g09_coords.py <file> (--input | --final) [--save <fname>] [--vmd]

Parse atomic coordinates from Gaussian input and output files

Argunemts:
    <file>             gjf file to parse
    --input            Print input xyz coords
    --final            Print last optimised xyz coords

Options:
    --save <fname>     Save coords to a file
    --vmd              Print in vmd style

pv278@cam.ac.uk, 07/08/15
"""
from docopt import docopt
import numpy as np
import re
from xyzlib import Atoms

def parse_input_coords(infile):
    """
    Parse the input atom names with xyz coords
    TODO: parse charge and multiplicity
    """
    regex = "Charge =(.*?)\\n \\n"
    with open(infile) as f:
        match = re.findall(regex, f.read(), re.S)   # WHAT IS re.S?
        match = match[-1].split("\n")[1:]
        match = [line.split() for line in match]
        names = [line[0] for line in match]
        coords = np.array([line[1:4] for line in match]).astype(float)
        return Atoms(names, coords)

def parse_last_coords(infile):
    """
    Parse last produced coords in a g09 simulation
    """
    names = parse_input_coords(infile).names
    regex = "Standard orientation(.*?)Rotational constants"
    with open(infile) as f:
        match = re.findall(regex, f.read(), re.S)
        match = match[-1].split("\n")[5:-2]
        match = [line.split() for line in match]
        coords = np.array(match).astype(float)[:, -3:]
        return Atoms(names, coords)


if __name__ == "__main__":
    args = docopt(__doc__)
#    print args
    infile = args["<file>"]

    if args["--input"]:
        A = parse_input_coords(infile)
        if args["--save"]:
            A.save(args["--save"], args["--vmd"])
        else:
            print A
    if args["--final"]:
        B = parse_last_coords(infile)
        if args["--save"]:
            B.save(args["--save"], args["--vmd"])
        else:
            print B


