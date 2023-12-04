#!/usr/bin/env python
"""Usage:
    default_ptfe_input_mix.py bulk [--seed <s> --L <L> --lmbda <l> --gamma <g>]
                   [--method <m> --r <r> --k0 <k0>]
    default_ptfe_input_mix.py slab [--width <w> --el <el> --rpt <rpt> --r <r>]
                   [--seed <s> --L <L> --lmbda <l> --gamma <g> --method <m> --k0 <k0>]
                          

Generate input.yaml file that will serve to generate
a LAMMPS/DL_MESO input files by being read by gen_ptfe.py script.
This script supersedes the default_ptfe_input.sh scripts.
New version with correct chi-params from solubilities only.

Options:
    --seed <s>     Random seed [default: 1234]
    --L <L>        Box size [default: 40]
    --lmbda <l>    Water uptake [default: 9]
    --gamma <g>    DPD friction [default: 4.5]
    --el <el>      Electrode type: carbon, silica, quartz [default: carbon]
    --width <w>    Nafion slab width in nm [default: 10]
    --rpt <rpt>    Ratio of width of Pt segment on electrode [default: 0]
    --method <m>   Use 25 (1), 6*25 (2) or 82.5 (3) as repulsion [default: 1]
    --r <r>        Polymerisation of Nafion chain [default: 15]
    --k0 <k0>      Bond coefficient [default: 4.0]

pv278@cam.ac.uk, 24/08/17
"""
import yaml
import sys
from docopt import docopt
from collections import OrderedDict


chi = {}
chi["carbon"] = \
"""
A B:   1.2300
A C:   7.4400
A E:   1.0958
A W:   3.3600
B C:   2.7000
B E:   0.9413
B W:   1.5300
C E:   0.0290
C W:   1.4800
E W:   3.7653
"""
chi["quartz"] = \
"""
A B:   1.2300
A C:   7.4400
A E:   3.6020
A W:   3.3600
B C:   2.7000
B E:   3.3171
B W:   1.5300
C E:   1.0430
C W:   1.4800
E W:   1.1867
"""


args = docopt(__doc__)
rc = 8.14e-10
seed = int(args["--seed"])
L = float(args["--L"])
T = 1.0
mono_beads = "AAABC"
Nmc = int(args["--r"])
lmbda = int(args["--lmbda"])
gamma = float(args["--gamma"])
k0 = float(args["--k0"])
r0 = 0.1
elmat = args["--el"].lower().rstrip("/")
w = float(args["--width"])
rPt = float(args["--rpt"])
method = int(args["--method"])

if elmat not in ["carbon", "silica", "quartz"]:
    sys.exit("ERROR. Choose from these electrodes: carbon, silica, quartz.")

L_nm = L * rc / 1e-9
if w > L_nm:
    sys.exit("Slab width can be at most the box size (%.2f nm)." % L_nm)

if rPt < 0.0 or rPt > 1.0:
    sys.exit("Platinum ratio must be between 0 and 1.")

s = """# ===== Bead types:
# * A, B, C
# * W: water bead: 6 H2O\n"""
s += "# * E: electrodes from %s\n" % elmat
s += """
# ===== Method:
# * 1: a_DPD = 25
# * 2: a_DPD = 158
# * 3: a_DPD = 82.5, from Fuchslin, JCP, 2009
\n"""

s += "seed:              %i\n" % seed
s += "method:            %i\n" % method

s += "box-size:          %.6f        # DPD units, 1 = 8.14 AA\n" % L
s += "temperature:       %.1f       # Units of kB T\n" % T
s += "mono-per-chain:    %i\n" % Nmc
s += "mono-beads:        %s\n" % mono_beads
s += "water-uptake:      %i         # number of H2O/SO3H\n" % lmbda
s += "gamma:             %.1f       # DPD drag coefficient\n\n" % gamma

s += "chi-params:\n"
chis = yaml.load(chi[elmat])
chis = OrderedDict(sorted(chis.items()))
#for k, v in yaml.load(chi[elmat]).items():
for k, v in chis.items():
    s += "    %s: %.3f\n" % (k, v)

s += "\n"
s += "bond-coeff:        %.1f\n" % k0
s += "equilibrium-dist:  %.1f\n\n" % r0

if args["slab"]:
    s += "electrodes:\n"
    s += "    material:      %s\n" % elmat
    s += "    width:         %.1f       # slab width in nm\n" % w 
    s += "    Pt-ratio:      %.1f       # ratio of Pt/C on electrode\n" % rPt

fname = "input.yaml"
open(fname, "w").write(s)
print("Parameter file saved in %s." % fname)


