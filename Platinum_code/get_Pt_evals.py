#!/usr/bin/env python
"""Usage:
    get_Pt_evals.py <cluster> [-d <d>] [--bandgap]

Arguments:
    <cluster>         Pt cluster, e.g. 9_10_9

Options:
    -h,--help         Show this message and exit
    -d <d>,--dir <d>  Directory of evals file
    --bandgap         Print band gap for a specific cluster for all spins

pv278@cam.ac.uk, 01/05/15
"""
from docopt import docopt
import numpy as np
import os
from itertools import chain
from iolib import save_table, print_table


def get_evals(cluster, Pt_dir, spin_list):
    Nall = 44*sum([int(i) for i in cluster.split("_")])
    Nocc = 18*sum([int(i) for i in cluster.split("_")])
    Nvirt = Nall - Nocc
    Ns = len(spin_list)
    mat_all = np.zeros((Nall, Ns))
    mat_occ = np.zeros((Nocc, Ns))
    mat_virt = np.zeros((Nvirt, Ns))

    for i in spin_list:
        outfile = open(get_path(Pt_dir, cluster, i)).readlines()
     
        all = [line.split()[2:] for line in [l for l in outfile if "Eigenvalues" in l]]
        e_all = [float(eval) for eval in list(chain(*all))]
     
        occ = [line.split()[4:] for line in [l for l in outfile if "occ." in l]]
        e_occ = [float(eval) for eval in list(chain(*occ))]
        
        virt = [line.split()[4:] for line in [l for l in outfile if "virt." in l]]
        e_virt = [float(eval) for eval in list(chain(*virt))]
     
        if i == 0:              # double num. of e-values due to degeneracy
            e_all += e_all
            e_occ += e_occ
            e_virt += e_virt
        
        print("S =",i,"\t",len(e_all),"\t",len(e_occ),"\t",len(e_virt))
    
        if len(e_all) == 0:     # if the run fails to converge
            e_all = [10] * Nall
            e_occ = [10] * Nocc
            e_virt = [10] * (Nall - Nocc)
    
        mat_all[:,i] = e_all
        mat_occ[:,i] = e_occ
        mat_virt[:,i] = e_virt
     
    return mat_all, mat_occ, mat_virt


def get_bandgap(cluster, spin, dirr="/home/pv278/Platinum/"):
    """Extract band gap (highest occ value - lowest virt value)
       and from the files for a specific spin"""
    outfile = open(get_path(dirr, cluster, spin)).readlines()
    line = [l for l in outfile if "occ" in l]
    if line:
        Eocc = line[-1].split()[-1]
    else:
        Eocc = None
    line = [l for l in outfile if "virt" in l]
    if line:
        Evirt = line[0].split()[4]
    else:
        Evirt = None
    
    if Eocc:
        return (float(Evirt) - float(Eocc)) * 27.211
    else:
        return


def get_all_bandgaps(cluster, spin_list, dir="/home/pv278/Platinum/"):
    """Extract band gaps for all spins"""
    E, s = [], []
    for spin in spin_list:
        Ebg = get_bandgap(cluster, spin, dir)
        if Ebg:
            E.append(Ebg)
            s.append(spin)
    return np.vstack((s,E)).T


if __name__ == "__main__":
    args = docopt(__doc__, version=1.0)
    
    cluster = args["<cluster>"]
    Pt_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..")) # "~/Platinum

    if args["--dir"]:
        base_dir = Pt_dir + "/" + args["--dir"]
    else:
        base_dir = Pt_dir + "/Plain"
    
    outdir = base_dir + "/Outfiles/Evals"
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        print(outdir, "created")
    
    spin_list = range(11)
    A_all, A_occ, A_virt = get_evals(cluster, Pt_dir, spin_list)
    
    outfile_all = outdir + "/evals_Pt" + cluster + ".out"
    outfile_occ = outdir + "/evals_Pt" + cluster + "_occ.out"
    outfile_virt = outdir + "/evals_pt" + cluster + "_virt.out"
    save_table(A_all, outfile_all)
    save_table(A_occ, outfile_occ)
    save_table(A_virt, outfile_virt)

    if args["--bandgap"]:
        A = get_all_bandgaps(cluster, spin_list)
        print_table(A)
        save_table(A,outdir + "/Pt" + cluster + "_bandgap.out")


