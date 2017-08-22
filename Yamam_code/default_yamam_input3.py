#!/usr/bin/env python
"""Usage:
    default_yamam_input.py bulk [--seed <s> --L <L> --lmbda <l> --gamma <g>]
    default_yamam_input.py slab [--width <w> --el <el> --r <r>]
               [--seed <s> --L <L> --lmbda <l> --gamma <g>]

Generate input.yaml file that will serve to generate
a LAMMPS data file by being read by gen_ptfe.py script.
This script supersedes the default_ptfe_input*.sh scripts.

Options:
    --seed <s>     Random seed [default: 1234]
    --L <L>        Box size [default: 40]
    --lmbda <l>    Water uptake [default: 9]
    --gamma <g>    DPD friction [default: 4.5]
    --el <el>      Electrode type: carbon, silica, quartz [default: carbon]
    --width <w>    Nafion slab width in nm [default: 10]
    --r <r>        Polymerisation of Nafion chain [default: 15]

pv278@cam.ac.uk, 05/07/17
"""
import yaml
import sys
from docopt import docopt


chi = {}
chi["carbon"] = \
"""
A B: 0.022
A C: 3.11
B C: 1.37
A W: 5.79
B W: 4.90
C W: -2.79
A E: 1.24
B E: 0.05
C E: 0.86
W E: 2.89
"""
chi["quartz"] = \
"""
A B: 0.022
A C: 3.11
B C: 1.37
A W: 5.79
B W: 4.90
C W: -2.79
A E: 1.82
B E: 0.25
C E: 0.25
W E: 0.94
"""


args = docopt(__doc__)
rc = 7.1e-10
seed = int(args["--seed"])
L = float(args["--L"])
T = 1.0
mono_beads = "AAAABC"
lmbda = float(args["--lmbda"])
Nmc = int(args["--r"])
gamma = 4.5
k0 = 4.0
r0 = 0.86
elmat = args["--el"].lower()
w = float(args["--width"])
rPt = 0.0
mono_beads = "AAAABC"

if elmat not in ["carbon", "silica", "quartz"]:
    print("ERROR. Choose from electrode support: carbon, silica, quartz.")
    sys.exit()

L_nm = L * rc / 1e-9
if w > L_nm:
    sys.exit("Slab width can be at most the box size (%.2f nm)." % L_nm)

if rPt < 0.0 or rPt > 1.0:
    sys.exit("ERROR: Platinum ratio must be between 0 and 1.")

s = """# ===== Bead types:
# * A, B, C
# * W: water bead: 6 H2O\n"""
s += "# * E: electrodes from %s\n" % elmat
s += """# * P: platinum
# =====\n"""

s += "seed:              %i\n" % seed
s += "box-size:          %.0f        # DPD units, 1 = 8.14 AA\n" % L
s += "temperature:       %.1f       # Units of kB T\n" % T
s += "mono-per-chain:    %i\n" % Nmc
s += "mono-beads:        %s\n" % mono_beads
s += "water-uptake:      %.2f       # \n" % lmbda
s += "gamma:             %.1f       # DPD drag coefficient\n" % gamma

s += "chi-params:\n"
for k, v in yaml.load(chi[elmat]).items():
    s += "    %s: %.2f\n" % (k, v)

s += "\n"
s += "bond-coeffs:\n    A A: %.1f\n\n" % k0
s += "equilibrium-dist:  %.2f\n\n" % r0

if args["slab"]:
    s += "electrodes:\n"
    s += "    material:      %s\n" % elmat
    s += "    width:         %.1f       # slab width in nm\n" % w 
    s += "    Pt-ratio:      %.1f       # ratio of Pt/C segment\n" % rPt

fname = "input.yaml"
open(fname, "w").write(s)
print("Parameter file saved in", fname)



