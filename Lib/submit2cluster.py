#!/usr/bin/env python
"""Usage:
    submit_cottrell.py gaussian <node> <dir> <fname> [--cores <nc> --dry --direct]
    submit_cottrell.py lammps <node> <dir> <fname> [--cores <nc> --dry]
    submit_cottrell.py dlms <node> [--cores <nc> --dry]
    submit_cottrell.py dlms_serial <node> [--dry]

Submit Gaussian, LAMMPS, or DL_MESO jobs to a cluster.
Calling submit2cottrell.sh with appropriate arguments.

Arguments:
    <node>        Cluster node, from 0 to 10
    <dir>         Directory of the gjf file
    <fname>       Name of the gjf file, NO EXTENSION
    --direct      Submit directly using Gaussian form

Options:
    --cores <nc>  Number of cores [default: 16]
    --dry         Dry run, print file cmd on screen

pv278@cam.ac.uk, 16/05/15, edit 08/08/16
"""
from docopt import docopt
import os, sys, subprocess


def get_sup_name(node, maxnum=27):
    """Get supervisor name"""
    if node > maxnum:
        sys.exit("Only nodes 0 to 27 available.")
    return "jae" if node < 10 else "pdb" if 10 <= node < 20 else "new"


if __name__ == "__main__":
    args = docopt(__doc__, version=0.1)
    node = int(args["<node>"])
    server = "%s.q@compute-0-%i" % (get_sup_name(node, 27), node)
#    server = get_sup_name(node, 27) + ".q@compute-0-" + str(node)
    cores = int(args["--cores"])

    ddir = args["<dir>"]
    filename  = args["<fname>"]
    bashscript = "/home/pv278/Res/PhDcode/Lib/submit2cluster.sh"

    if args["gaussian"]:
        prog = "gaussian"
        filepath = os.path.join(os.getcwd(), ddir, filename)
        filepath = filepath.rstrip(".gjf")
        if args["--direct"]:
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
            submit_string = "qsub -b y -q %s -pe orte %i %s" % \
                (server, cores, cmd)
#            submit_string = "qsub -b y -q " + server + " -pe orte " + str(cores) + " " + cmd
 
        else:
            submit_string = "qsub -q %s -pe orte %i %s %s %s" % \
                (server, cores, bashscript, prog, filepath)
#            submit_string = "qsub -q " + server + " -pe orte " + str(cores) + \
#                            " " + bashscript + " " + prog + " " + filepath

    elif args["lammps"]:
        prog = "lammps"
        filepath = os.path.join(os.getcwd(), ddir, filename)
        submit_string = "qsub -q %s -pe orte %i %s %s %s %i" % \
            (server, cores, bashscript, prog, filepath, cores)

    elif args["dlms"]:
        prog = "dlms"
        filepath = "dummy"
        submit_string = "qsub -q %s -pe orte %i %s %s %s %i" % \
            (server, cores, bashscript, prog, filepath, cores)

    elif args["dlms_serial"]:
        prog = "dlms_serial"
        filepath = "dummy"
        cores = 1
        submit_string = "qsub -q %s -pe orte %i %s %s %s %i" % \
            (server, cores, bashscript, prog, filepath, cores)


    if args["--dry"]:
        print(submit_string)
    else:
        subprocess.call(submit_string, shell=True)

