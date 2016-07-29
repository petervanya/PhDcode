#!/usr/bin/python
"""
Plotting density of states for Pt atoms in one layer
Created: 21/01/14
"""
import matplotlib.pyplot as plt
import numpy as np
import argparse

# ===== parse arguments
help="Plot density of state for a given cluster."
parser=argparse.ArgumentParser(description=help,epilog="Author: pv278@cam.ac.uk")

parser.add_argument("-fi","--fin",dest="fin",action="store",type=str,required=True,
                    metavar="f",help="Input file for e-values")

parser.add_argument("-fioc","--fin_occ",dest="fin_occ",action="store",type=str,required=True,
                    metavar="f",help="Input file for occupied e-values")

parser.add_argument("-fo","--fout",dest="fout",action="store",type=str,
                    metavar="f",help="Output file (not required)")
                                        
parser.add_argument("-s","--spin",dest="spin",action="store",type=int,default=0,
                    metavar="s",help="Spin of the state to plot")

args=parser.parse_args()

# ===== plotting
E = np.loadtxt(args.fin)
Eocc = np.loadtxt(args.fin_occ)
spin=args.spin
print "Spin =",spin

#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')

xtics=[-4,-3,-2,-1,0,1]
xlbl="$E-E_\mathrm{F}$ (Hartree)"
ylbl="DoS"
xmin=-4
xmax=1

Bins=np.arange(xmin,xmax+0.2,0.1)
plt.hist(E[:,spin],    bins=Bins,edgecolor="blue")
plt.hist(Eocc[:,spin], bins=Bins,color="red",edgecolor="red")

plt.xlabel(xlbl)
plt.ylabel(ylbl)
plt.xticks(xtics)
plt.xlim([xmin,xmax])
plt.ylim([0,20])
plt.title("Density of states")

#plt.tight_layout()
# ===== save figure
if args.fout:
  outfile=args.fout
else:
  outfile=args.fin.replace(args.fin.split(".")[-1],"pdf")
  
plt.savefig(outfile)
print "Figure saved in",outfile
