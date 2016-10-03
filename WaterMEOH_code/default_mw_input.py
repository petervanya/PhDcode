#!/usr/bin/env python
"""Usage:
    default_mw_input.py [options]
    
[--rho <rho> --L <L> --T <T> --dt <dt>]
[--runtime <rt> --fw <fw> --method <m> --bead <bw> <bm>]

Options:
    --rho <rho>         DPD density [default: 3]
    --L <L>             Box size in nm [default: 10]
    --T <T>             Temperature in K [default: 300]
    --dt <dt>           Time step in DPD units [default: 0.05]
    --runtime <rt>      Total run time in ns [default: 10]
    --fw <fw>           Molar fraction of water [default: 0.5]
    --method <m>        CG method [default: 1]
    --bead-mols <bm>    Bead compositions, e.g. "3 4" [default: 2 1]

pv278@cam.ac.uk, 08/08/16
"""
from docopt import docopt


args = docopt(__doc__)
rho = float(args["--rho"])
L = float(args["--L"])
T = float(args["--T"])
dt = float(args["--dt"])
runtime = float(args["--runtime"])
fw = float(args["--fw"])
method = int(args["--method"])
bw, bm = map(int, args["--bead-mols"].split())


s = """# Methods:
# 1. aii = 25 kT, da ~chi from FH
# 2. aii based on compressibility, da = sqrt(a11 a22) + ~chi
# 3. from Travis, JCP, 2007
"""

s += "dpd-density: %.1f\n" % rho
s += "box-size: %.1f             # in nm\n" % L
s += "temperature: %.1f         # in K\n" % T
s += "dt: %.3f                   # in DPD units\n" % dt
s += "run-time: %.1f             # in ns\n" % runtime
s += "water-molar-fraction: %.2f # fraction of water molecules\n" % fw
s += "method: %i\n" % method
s += "waters-in-bead: %i\n" % bw
s += "meohs-in-bead: %i\n" % bm

fname = "input.yaml"
open(fname, "w").write(s)
print("Parameter file saved in %s." % fname)

