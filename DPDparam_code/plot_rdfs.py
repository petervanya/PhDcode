#!/usr/bin/env python
"""Usage:
    plot_rdfs.py <rdffiles> [--title <t>]

Plot all RDFs in one figure.

Options:
    --title <t>    Title [default: "title"]

05/04/17
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import glob
from docopt import docopt

matplotlib.rcParams.update({"font.size": 12})

def parse_nm(s):
    tmp = s.split("/")[-1].rstrip(".out").split("_")
    hit = [elem for elem in tmp if "nm" in elem]
    return int(hit[-1].lstrip("nm"))


args = docopt(__doc__)
files = glob.glob(args["<rdffiles>"])
files.sort(key=parse_nm)

for f in files:
    A = np.loadtxt(f)
    nm = parse_nm(f)
    lbl = "$N_{\\mathrm{m}}=%i$" % nm
    plt.plot(A[:, 0], A[:, 1], lw=2, label=lbl)

plt.xlim([0, 20])
plt.yticks(np.linspace(0, 1.5, 6))
plt.ylim([0, 1.5])
plt.xlabel("$r$ (A)")
plt.ylabel("$g(r)$")
plt.legend(loc=4, ncol=2)
plt.title(args["--title"].rstrip("/"))

figname = "rdf_all.png"
plt.savefig(figname, bbox_inches="tight")
print("Master figure saved in %s." % figname)

