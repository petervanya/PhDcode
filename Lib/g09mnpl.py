#!/usr/bin/env python
"""Usage:
    g09mnpl.py gen_header [--nproc <np>] [--method <m>] [--basis <bs>] 
                          [--opt] [--freq] [--scfsteps <scf>] [--misc <misc>]
    g09mnpl.py parse_coords <file> (--input | --final)
    g09mnpl.py gen_gjf from_xyz <files> [--header <h>] [--save]
    g09mnpl.py gen_gjf from_out <file> <newname>
    g09mnpl.py extract energies <files>  [--save <datafile>]

Gaussian 09 manipulate.
Collection of method to manipulate files entering or leaving Gaussian.
* construct gjf files
* parse xyz coordinates from outfiles
* generate headers
* extract data from outfiles

Options:
    --nproc <np>         Number of processors to use [default: 16]
    --method <m>         Qchem method [default: B3LYP]
    --basis <bs>         Basis set [default: 6-31g*]
    --opt                Produce "Opt" flag
    --freq               Produce "Freq" flag
    --scfsteps <scf>     Number of self-consistent steps [default: 1000]
    --misc <misc>        Various other keywords
    --header <h>         Set Gaussian header

pv278@cam.ac.uk, 06/10/15
"""
import numpy as np
import pandas as pd
import glob, os, sys, re
from docopt import docopt
from xyz_lib import Atoms


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
    """Parse last produced coords in a g09 simulation"""
    names = parse_input_coords(infile).names
    regex = "Standard orientation(.*?)Rotational constants"
    with open(infile) as f:
        match = re.findall(regex, f.read(), re.S)
        match = match[-1].split("\n")[5:-2]
        match = [line.split() for line in match]
        coords = np.array(match).astype(float)[:, -3:]
        return Atoms(names, coords)


def gen_header(p):
    """Generate Gaussian header, arguments: method, basis, opt"""
    header = "%nproc=" + str(p.get("np")) + "\n"
    header += "#T " + p.get("method") + "/" + p.get("basis") \
           + " Test " + p.get("opt") + " " + p.get("freq") \
           + " " + p.get("scf") + " " + p.get("misc")
    header += "\n\n" + "Blabla" + "\n\n"
    return header

 
def get_params(args):
    """Get parameters from command line arguments
    for Gaussian header"""
    params = {}
    params["np"] = args["--nproc"]
    params["method"] = args["--method"]
    params["basis"] = args["--basis"]
    params["opt"] = "Opt" if args["--opt"] else ""
    params["freq"] = "Freq" if args["--freq"] else ""
    params["scf"] = "scf=(direct, maxcycle=%s)" % args["--scfsteps"] \
            if args["--scfsteps"] else ""
    params["misc"] = args["--misc"] if args["--misc"] else ""
    return params 


if __name__ == "__main__":
    args = docopt(__doc__)
    default_header = "%nproc=16\n#T B3LYP/6-31G* Test\n\nBlabla\n\n"  # AD HOC SOL
    params = get_params(args)

    if args["gen_header"]:
        header = gen_header(params)
        print(header)
    
    if args["parse_coords"]:    #TO TEST
        infile = args["<file>"]
        if args["--input"]:
            A = parse_input_coords(infile)
            if args["--save"]:
                A.save(args["--save"], args["--vmd"])
            else:
                print(A)
        if args["--final"]:
            B = parse_last_coords(infile)
            if args["--save"]:
                B.save(args["--save"], args["--vmd"])
            else:
                print(B)

    elif args["gen_gjf"]:
        if args["from_xyz"]:
            xyzfile = args["<files>"]
            A = Atoms().read(xyzfile)
            if args["--header"]:
                string = args["--header"] + "\n\n" + str(A) + "\n\n"
            else:
                string = default_header + str(A) + "\n\n"
            fname = xyzfile.rstrip("xyz") + "gjf"
            open(fname, "w").write(string)
            print("gjf file written into %s." % fname)

#        if args["from_xyz"]:   # MAKE THIS QUICKER, NOW TOO SLOW
#            xyzfiles = glob.glob(args["<files>"])
##            xyzfiles = args["<files>"]
#            A = Atoms().read(xyzfiles[0])
#            if len(xyzfiles) > 0:
#                for xyzfile in xyzfiles[1:]:
#                    A = A + Atoms().read(xyzfile)
#            if args["--header"]:
#                string = args["--header"] + "\n\n" + str(A) + "\n\n"
#            else:
#                string = default_header + str(A) + "\n\n"
#            fname = xyzfiles[0].rstrip("xyz") + "gjf"
#            open(fname, "w").write(string)
#            print("gjf file written into", fname)

        if args["from_out"]:  #TO TEST, EXTRACT HEADER FROM OUTFILE
            outfile = args["<file>"]
            assert(outfile[-3:] == "out")   # NOT USER FRIENDLY
            A = parse_last_coords(outfile)
            text = open(outfile, "r").readlines()
            nproc = [line for line in text if "%nproc" in line][0][1:]  # skipping first blank character
            route = [line for line in text if " #" in line][0][1:]      # same
            #print re.findall(r"(^.*?%s.*?$)" % "nproc", open(outfile).read(), re.MULTILINE)
            # print re.findall(r"(^.*?%s.*?$)" % " #", open(outfile).read(), re.MULTILINE)
            header = nproc + route + "\nBlabla\n\n"
            newfname = args["<newname>"]
            open(newfname, "w").write(header + str(A) + "\n\n")
            print("gjf file written into", newfname)

    if args["extract"] and args["energies"]:
        outfiles = glob.glob(args["<files>"])
        print(outfiles)
        positions, energies = [], []
        
        for outfile in outfiles:
            outfile = os.getcwd() + "/" + outfile
            E_line = [line for line in open(outfile, "r").readlines() \
                    if "SCF Done" in line]
            if E_line:       # if energy line exists
                E_line = E_line[0]
                energies.append(float(E_line.split()[4]))
                # find distance in file name, NOT GENERAL FOR THIS SCRIPT
                pos = re.search(r"(?<=_d)[^}]*(?=.out)", outfile)
                pos = pos.group()
                positions.append(float(pos))
        
        if not energies:     # if the array is empty
            sys.exit()
        tbl = pd.DataFrame(energies, index=positions).sort_index()
        tbl.columns = ["E"]
        print(tbl)
        if args["<datafile>"]:
            datapath = args["<datafile>"]
            tbl.to_csv(datapath, sep="\t", header=False)
            print("Table saved in ", datapath)


