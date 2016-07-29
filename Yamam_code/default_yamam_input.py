#!/usr/bin/env python
"""Usage:
    default_ptfe_input.py [--L <L> --wt <wt>]
                          [--el <el> --width <w>]

Generate input.yaml file that will serve to generate
a LAMMPS data file by being read by gen_ptfe.py script.
This script supersedes the default_ptfe_input*.sh scripts.

Options:
    --L <L>             Box size [default: 40]
    --wt <wt>           Water content [default: 0.2]

pv278@cam.ac.uk, 16/06/16
"""
import yaml
import sys
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
L = float(args["--L"])
T = 1.0
Nmc = 5            # very different!
wt = float(args["--wt"])
gamma = 4.5
k0 = 100.0
r0 = 0.86
elmat = "carbon"
w = 0.0
rPt = 0.0

if elmat not in ["carbon", "silica", "quartz"]:
    print("ERROR. Choose from electrode support: carbon, silica, quartz.")
    sys.exit()

if wt < 0.0 or wt > 1.0:
    sys.exit("Water content must be between 0 and 1.")

if w < 0.0 or w > L/2:
    print("ERROR: Electrode width is at most one half of box size.")
    sys.exit()

if rPt < 0.0 or rPt > 1.0:
    print("ERROR: Platinum ratio must be between 0 and 1.")
    sys.exit()

s = """# ===== Bead types:
# * A, B, C
# * W: water bead: 6 H2O\n"""
s += "# * E: electrodes from %s\n" % elmat
s += """# * P: platinum
# =====\n"""

s += "box-size:          %.0f        # DPD units, 1 = 8.14 AA\n" % L
s += "temperature:       %.1f       # Units of kB T\n" % T
s += "mono-per-chain:    %i\n" % Nmc
s += "water-content:     %.2f       # ratio < 1\n" % wt
s += "gamma:             %.1f       # DPD drag coefficient\n" % gamma
s += "topology:          (A 3, [B 1], [C 1])15\n\n"

s += "chi-params:\n"
for k, v in yaml.load(chi[elmat]).items():
    s += "    %s: %.2f\n" % (k, v)

s += "\n"
s += "bond-coeffs:\n    A A: %.1f\n\n" % k0
s += "equilibrium-dist:  %.1f\n\n" % r0

s += "electrodes:\n"
s += "    material:      %s\n" % elmat
s += "    width:         %.1f       # electrode width\n" % w 
s += "    Pt-ratio:      %.1f       # ratio of Pt/C segment\n" % rPt

fname = "input.yaml"
open(fname, "w").write(s)
print("Parameter file saved in", fname)



