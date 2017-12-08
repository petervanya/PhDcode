#!/usr/bin/env python
"""Usage:
    default_mdpd_ptfe.py [--we <we> --wf <wf> --el <el> --rpt <rpt>]
                         [--seed <s> --L <L> --lmbda <l> --gamma <g>]
                         [--A <A> --B <B>]
                          

Create input.yaml file that will serve to generate
a LAMMPS/DL_MESO input files by being read by gen_ptfe.py script.
This script supersedes the default_ptfe_input*.sh scripts.
Default A and B derived for Nm = 6 by fitting density and diffusivity.

Options:
    --seed <s>     Random seed [default: 1234]
    --L <L>        Box size [default: 40]
    --lmbda <l>    Water uptake [default: 9]
    --gamma <g>    DPD friction [default: 15]
    --el <el>      Electrode type: carbon, silica, quartz [default: carbon]
    --we <we>      Electrode width in nm [default: 10]
    --wf <wf>      Nafion film width in nm [default: 10]
    --rpt <rpt>    Ratio of width of Pt segment on electrode [default: 0]
    --A <A>        MDPD repulsion [default: -80.4]
    --B <B>        MDPD attraction [default: 51.7]

pv278@cam.ac.uk, 04/12/17
"""
import yaml
import sys
from docopt import docopt


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
Nmc = 15
lmbda = int(args["--lmbda"])
gamma = float(args["--gamma"])
k0 = 4.0
r0 = 0.1
elmat = args["--el"].lower()
we = float(args["--we"])
wf = float(args["--wf"])
rPt = float(args["--rpt"])
A_def = float(args["--A"])
B_def = float(args["--B"])
rho = 2.901 + 0.68 * (-A_def) * (B_def - 4.09)**(-0.716)
chi2da = - 0.259 + 0.196 * rho 

if elmat not in ["carbon", "quartz"]:
    sys.exit("ERROR. Choose from these electrodes: carbon, quartz.")

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
s += "# * P: platinum\n"

s += "seed:              %i\n" % seed
s += "A:                 %.2f       # interaction param\n" % A_def
s += "B:                 %.2f       # interaction param\n" % B_def
s += "density:           %.3f       # calculated from A, B\n" % rho
s += "chi2da:            %.3f       # linear fit of Jamali, JCP, 2015\n" % chi2da

s += "box-size:          %.0f       # DPD units\n" % L
s += "temperature:       %.1f       # Units of kT\n" % T
s += "polymerisation:    %i\n" % Nmc
s += "mono-beads:        %s\n" % mono_beads
s += "water-uptake:      %i         # number of H2O/SO3H\n" % lmbda
s += "gamma:             %.1f       # DPD friction\n" % gamma

s += "chi-params:\n"
for k, v in yaml.load(chi[elmat]).items():
    s += "    %s: %.2f\n" % (k, v)

s += "\n"
s += "bond-coeff:        %.1f\n" % k0
s += "equilibrium-dist:  %.1f\n\n" % r0

s += "film-width:    %.2f       # film width in nm\n" % wf
s += "el-mat:        %s\n" % elmat
s += "el-width:      %.2f       # electrode width in nm\n" % we
s += "Pt-ratio:      %.1f       # ratio of Pt/C on electrode\n" % rPt

fname = "input.yaml"
open(fname, "w").write(s)
print("Parameter file saved in %s." % fname)


