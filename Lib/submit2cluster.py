#!/usr/bin/env python
"""Usage:
    submit_cottrell.py gaussian <node> <dir> <fname> (--bash | --direct) [--cores <nc>] [--dry]
    submit_cottrell.py lammps <node> <dir> <fname> (--bash | --direct) [--cores <nc>] [--dry]

A script to easily submit Gaussian and LAMMPS jobs to the Cottrell cluster

Arguments:
    <node>        Cluster node, from 0 to 10
    <dir>         Directory of the gjf file
    <fname>       Name of the gjf file, NO EXTENSION
    --direct      Submit directly using Gaussian form
    --bash        Submit via bash script cottrell.sh

Options:
    --cores <nc>  Number of cores [default: 16]
    --dry         Dry run, print file cmd on screen

pv278@cam.ac.uk, 16/05/15
"""
from docopt import docopt
import os, sys, subprocess

def get_sup(node, maxnum=27):
    """Get supervisor name"""
    if node > maxnum:
        print "Only nodes 0 to 27 available, aborting."
        sys.exit()
    return "jae" if node < 10 else "pdb" if 10 <= node < 20 else "new"

if __name__ == "__main__":
    args = docopt(__doc__, version=0.1)
#    print args
    node = int(args["<node>"])
    server = get_sup(node, 27) + ".q@compute-0-" + str(node)
    cores = int(args["--cores"])

    ddir = args["<dir>"]
    filename  = args["<fname>"]
    bashscript = "/home/pv278/GeneralScripts/submit2cluster.sh"

    if args["gaussian"]:
        prog = "gaussian"
        filepath = os.path.join(os.getcwd(), ddir, filename)
        filepath = filepath.rstrip(".gjf")
        if args["--bash"]:
            submit_string = "qsub -q " + server + " -pe orte " + str(cores) + \
                            " " + bashscript + " " + prog + " " + filepath
            if args["--dry"]:
                print submit_string
                sys.exit()
            subprocess.call(submit_string, shell=True)
        elif args["--direct"]:
            ext = ".gjf"
            base_dir = os.path.abspath(os.path.join(os.getcwd(), ".."))
            filedir = os.path.join(base_dir, dir)
            infilepath = filedir + "/" + filename + ext
            outfilepath = filedir + "/" + filename + ".out"
            
            exports = """
            g09root=\"/home/Gaussian/\"
            GAUSS_SCRDIR=\"/state/partition1/Gaussian_scratch\"
            GAUSS_EXEDIR=\"/home/Gaussian/g09/bsd:/home/Gaussian/g09/private:/home/Gaussian/g09\"
            export g09root GAUSS_SCRDIR GAUSS_EXEDIR
            . $g09root/g09/bsd/g09.profile"""
            subprocess.call(exports, shell=True)
            
            gaussianbin = "/home/Gaussian/g09/g09"
            cmd = gaussianbin + " < " + infilepath + " > " + outfilepath
            submit_string = "qsub -b y -q " + server + " -pe orte " + str(cores) + " " + cmd
 
            if args["--dry"]:
                print submit_string
                sys.exit()
            subprocess.call(submit_string, shell=True)

    elif args["lammps"]:
        prog = "lammps"
        filepath = os.path.join(os.getcwd(), ddir, filename)
        submit_string = "qsub -q " + server + " -pe orte " + str(cores) + \
                         " " + bashscript + " " + \
                         prog + " " + filepath + " " + str(cores)
        if args["--dry"]:
            print submit_string
            sys.exit()
        subprocess.call(submit_string, shell=True)


