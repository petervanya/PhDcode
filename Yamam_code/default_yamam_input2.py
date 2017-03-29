#!/usr/bin/env python
"""Usage:
    default_yamam_input2.py [--L <L> --lmbda <l> --seed <s> --r <r>]

Generate input.yaml file that will serve to generate
a LAMMPS data file by being read by gen_ptfe.py script.
This script supersedes the default_ptfe_input*.sh scripts.

Options:
    --seed <s>     Random seed [default: 1]
    --lmbda <l>    Water uptake [default: 9]
    --L <L>        Box size [default: 40]
    --r <r>        Polymerisation [default: 5]

pv278@cam.ac.uk, 22/02/17
"""
import yaml
from docopt import docopt


chi = {}
chi["carbon"] = \
"""
# adsorption.xlsx
A B: 0.022
A C: 3.11
B C: 1.37
A W: 5.79
B W: 4.90
C W: -2.79
"""


args = docopt(__doc__)
seed = int(args["--seed"])
L = float(args["--L"])
T = 1.0
Nmc = int(args["--r"])         # very different!
lmbda = float(args["--lmbda"])
gamma = 4.5
k0 = 100.0
r0 = 0.86
elmat = "carbon"
w = 0.0
rPt = 0.0

s = """# ===== Bead types:
# * A, B, C
# * W: water bead: 6 H2O\n"""
s += "# * E: electrodes from %s\n" % elmat
s += """# * P: platinum
# =====\n"""

s += "seed:              %i\n" % seed
s += "box-size:          %.0f        # DPD units\n" % L
s += "temperature:       %.1f       # Units of kB T\n" % T
s += "mono-per-chain:    %i\n" % Nmc
s += "water-uptake:      %i\n" % lmbda
s += "gamma:             %.1f       # DPD drag coefficient\n" % gamma

s += "chi-params:\n"
for k, v in yaml.load(chi[elmat]).items():
    s += "    %s: %.2f\n" % (k, v)

s += "\n"
s += "bond-coeffs:\n    A A: %.1f\n\n" % k0
s += "equilibrium-dist:  %.2f\n\n" % r0

s += "electrodes:\n"
s += "    material:      %s\n" % elmat
s += "    width:         %.1f       # electrode width\n" % w 
s += "    Pt-ratio:      %.1f       # ratio of Pt/C segment\n" % rPt

fname = "input.yaml"
open(fname, "w").write(s)
print("Parameter file saved in", fname)



