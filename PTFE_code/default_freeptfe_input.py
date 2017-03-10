#!/usr/bin/env python
"""Usage:
    default_ptfe_input.py [--we <we> --wf <wf> --el <el> --rpt <rpt>]
                   [--seed <s> --L <L> --lmbda <l> --gamma <g> --method <m>]
                          

Generate input.yaml file that will serve to generate
a LAMMPS/DL_MESO input files by being read by gen_ptfe.py script.
This script supersedes the default_ptfe_input*.sh scripts.

Options:
    --seed <s>     Random seed [default: 1234]
    --L <L>        Box size [default: 40]
    --lmbda <l>    Water uptake [default: 9]
    --gamma <g>    DPD friction [default: 4.5]
    --el <el>      Electrode type: carbon, silica, quartz [default: carbon]
    --we <we>      Electrode width in nm [default: 10]
    --wf <wf>      Nafion film width in nm [default: 10]
    --rpt <rpt>    Ratio of width of Pt segment on electrode [default: 0]
    --method <m>   Use either 25 (1) or 6*25 (2) as repulsion [default: 1]

pv278@cam.ac.uk, 31/01/16
"""
import yaml
import sys
from docopt import docopt

"""
# My calculation for carbon
# use calculation below by JAE from adsorption.xlsx
A B: 1.23
A C: 7.44
B C: 2.70
A W: 3.36
B W: 1.53
C W: 1.48
A E: 0.29
B E: 0.06
C E: 0.06
W E: 2.96
"""

chi = {}
chi["carbon"] = \
"""
# adsorption.xlsx
A B: 1.23
A C: 7.44
B C: 2.53
A W: 16.91
B W: 8.66
C W: 1.85
A E: 1.24
B E: 0.05
C E: 0.86
W E: 2.89
A P: 15.77
B P: 8.07
C P: 1.82
W P: 0.01
E P: 1.98
"""
chi["quartz"] = \
"""
A B: 1.23
A C: 7.44
B C: 2.70
A W: 3.36
B W: 1.53
C W: 1.48
A E: 1.82
B E: 0.25
C E: 0.25
W E: 0.94
"""
chi["silica"] = \
"""
A B: 1.23
A C: 7.44
B C: 2.70
A W: 3.36
B W: 1.53
C W: 1.48
A E: 1.51
B E: 0.10
C E: 0.10
W E: 1.56
"""


args = docopt(__doc__)
rc = 8.14e-10
seed = int(args["--seed"])
L = float(args["--L"])
T = 1.0
mono_beads = "AAABC"
Nmc = 15
lmbda = int(args["--lmbda"])
gamma = float(args["--gamma"])
k0 = 4.0
r0 = 0.1
elmat = args["--el"].lower()
we = float(args["--we"])
wf = float(args["--wf"])
rPt = float(args["--rpt"])
method = int(args["--method"])

if elmat not in ["carbon", "silica", "quartz"]:
    sys.exit("ERROR. Choose from these electrodes: carbon, silica, quartz.")

#if w < 0.0 or w > L/2:
#    sys.exit("ERROR: Electrode width is at most one half of box size.")

L_nm = L * rc / 1e-9
if (we + wf) > L_nm:
    sys.exit("Slab width can be at most the box size (%.2f nm)." % L_nm)

if rPt < 0.0 or rPt > 1.0:
    sys.exit("Platinum ratio must be between 0 and 1.")

s = """# Input file for free space simulations.
# ===== Bead types:
# * A, B, C
# * W: water bead: 6 H2O\n"""
s += "# * E: electrodes from %s\n" % elmat
s += "# * P: platinum"
s += """
# ===== Method:
# * 1: a_DPD = 25 = (16 - 1) / 0.2 / 3
# * 2: a_DPD = 158 = (6*16 - 1) / 0.2 / 3
# * 3: a_DPD sim N^(2/3), from Fuchslin, JCP, 2009
\n"""

s += "seed:              %i\n" % seed
s += "method:            %i\n" % method

s += "box-size:          %.0f        # DPD units, 1 = 8.14 AA\n" % L
s += "temperature:       %.1f       # Units of kB T\n" % T
s += "polymerisation:    %i\n" % Nmc
s += "mono-beads:        %s\n" % mono_beads
s += "water-uptake:      %i         # number of H2O/SO3H\n" % lmbda
s += "gamma:             %.1f       # DPD drag coefficient\n" % gamma

s += "chi-params:\n"
for k, v in yaml.load(chi[elmat]).items():
    s += "    %s: %.2f\n" % (k, v)

s += "\n"
s += "bond-coeffs:\n    A A: %.1f\n\n" % k0
s += "equilibrium-dist:  %.1f\n\n" % r0

s += "film-width:    %.1f       # film width in nm\n" % wf
s += "el-mat:        %s\n" % elmat
s += "el-width:      %.1f       # electrode width in nm\n" % we
s += "Pt-ratio:      %.1f       # ratio of Pt/C on electrode\n" % rPt

fname = "input.yaml"
open(fname, "w").write(s)
print("Parameter file saved in %s." % fname)


