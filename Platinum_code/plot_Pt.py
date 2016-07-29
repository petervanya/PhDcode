#!/usr/bin/env python
"""Usage:
    plot.py dos <cluster> [-s <s>] [--convolved [--sigma <s>]]
    plot.py bandgap <cluster>

Plot density of states or bandgap for Pt atoms and save figures
in Outfiles/Evals

Arguments:
    <cluster>           Pt cluster, e.g. 9_10_9

Options:
    -s <s>,--spin <s>   Spin state to plot
    --convolved         Continuous DoS plot using gaussians instead of histograms
    --sigma <s>         Std dev of the convoluting gaussian [default: 1.0]

pv278@cam.ac.uk, 21/01/14
"""
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from math import exp, sqrt, pi
import re
from docopt import docopt

def plot_spin(E1, E2, fout):
    #plt.rc('text', usetex=True)
    #plt.rc('font', family='serif')
    xmin = -10
    xmax = 30
    
    Bins = np.arange(xmin, xmax+1, 0.3)
    plt.hist(E1, bins=Bins, edgecolor="blue")
    plt.hist(E2, bins=Bins, color="red", edgecolor="red")
    
    plt.xlabel("$E - E_F$ (eV)")
    plt.ylabel("DoS")
    plt.xlim([xmin, xmax])
    spin = re.findall("S\d{1,2}", fout)[0][1:]
    plt.title("Density of states, $S$=" + spin)

    plt.savefig(fout)
    plt.close()
    print "Figure saved in", fout

def plot_spin_convolved(MO, occMO, fout, sigma):
    """Plot DoS with gaussian-smoothed MOs"""
    fs = 16
    xmin = min(MO) - 1.0
    xmax = max(MO) + 1.0
    xlim = [-10, 30]
    delta = 0.1
    energy = np.arange(xmin, xmax, delta)
    DoS = [0]*len(energy)
    for E in MO:
        mid = np.where(energy == min(energy, key=lambda x: abs(x-E)))[0][0]
        for i in range(mid-10, mid+10):
            DoS[i] += gaussian(energy[i], E, sigma)

    #fig, axes = plt.subplots(2, 1, sharex=True)
    fig = plt.figure()
    
    ax1 = fig.add_subplot(2,1,1)
    ax1.plot(energy, DoS, linewidth=2.0)
    #ax1.set_ylabel("DoS", fontsize=fs)
    ax1.set_xlim(xlim)
    spin = re.findall("S\d{1,2}", fout)[0][1:]
    ax1.set_title("Density of states, $S$=" + spin, fontsize=fs)
    ax1.set_xticklabels([])
    ax1.set_yticklabels([])

    ax2 = fig.add_subplot(2,1,2)
    Bins = np.arange(xlim[0], xlim[1], 0.01)
    ax2.hist(MO, bins=Bins, color="r", edgecolor="r")
    ax2.hist(occMO, bins=Bins, color="g", edgecolor="g")
    ax2.set_xlabel("$E$ - $E_F$ (eV)", fontsize=fs)
    ax2.set_ylim([0, 1])
    ax2.set_yticks([0,1])
    ax2.set_yticklabels([])
    ax2.set_aspect(4.0)

    fig.subplots_adjust(hspace=-0.4)

    plt.savefig(fout)
    plt.close()
    print "Figure saved in", fout


def plot_bandgap(cluster, dir):
    """function to plot bandgap v.s. spin state"""
    fs = 16
    fin = dir + "Pt" + cluster + "_bandgap.out"
    A = np.loadtxt(fin)
    spin = A[:,0]
    Eg = A[:,1]

    sct = plt.plot(spin, Eg, "bo", ms=10)
#    plt.setp(sct, lw=4.0)
    plt.xlim([-0.1, 10.1])
    plt.xlabel("Spin", fontsize=fs)
    plt.ylabel("$E_{\\mathrm{gap}}$ (eV)", fontsize=fs)
    plt.title("Pt " + cluster.replace("_","."), fontsize=fs)

    fout = fin.replace(".out", ".png")
    plt.savefig(fout)
    plt.close()
    print "Figure saved in",fout


def gaussian(x, mu, sigma=0.1):
    return exp(-(x-mu)**2/(2*sigma**2)) / sqrt(2*pi*sigma**2)


def shift_to_Fermi(E, Eocc):
    """Shift e-values so that E=0 at the Fermi surface"""
    HOMO = Eocc[-1]
    E -= HOMO
    Eocc -= HOMO
    return E, Eocc


if __name__ == "__main__":
    args = docopt(__doc__, version=1.0)
#    print args
    cluster = args["<cluster>"]
    base_dir = "/home/pv278/Platinum/Plain/Outfiles/Evals/"
    
    if args["dos"]:
        fin = base_dir + "evals_Pt" + cluster + ".out"
        fin_occ = base_dir + "evals_Pt" + cluster + "_occ.out"
        spin = args["--spin"]
    
        E = np.loadtxt(fin)
        Eocc = np.loadtxt(fin_occ)
        Ha = 27.211
        E, Eocc = E*Ha, Eocc*Ha      # conversion to eV
        E, Eocc = shift_to_Fermi(E, Eocc)
        
        if spin:
            spin_list = [int(spin)]
        else:
            spin_list = range(11)
        
        for spin in spin_list:
            fout = fin.replace(".out", "_S" + str(spin) + ".png")
            if args["--convolved"]:
                sigma = float(args["--sigma"])
                plot_spin_convolved(E[:,spin], Eocc[:,spin], fout, sigma)
            else:
                plot_spin(E[:,spin], Eocc[:,spin], fout)

    if args["bandgap"]:
        plot_bandgap(cluster, base_dir)


