#!/usr/bin/env python
"""Usage:
    dlms_control.py [--L <L> --dt <dt> --steps <n> --thermo <th> --halo <h>]

Generate DL_MESO control file.

Options:
     --L <L>            Box length [default: 40.0]
     --dt <dt>          Timestep [default: 0.05]
     --steps <n>        Number of steps [default: 10000]
     --thermo <th>      Print every [default: 100]
     --halo <h>         Boundary halo, like neighbor [default: 2.5]

pv278@cam.ac.uk, 06/06/16
"""
from docopt import docopt
import sys


args = docopt(__doc__)
L = float(args["--L"])
dt = float(args["--dt"])
N = int(args["--steps"])
thermo = int(args["--thermo"])
halo = float(args["--halo"])


s = "pokus\n\n"

s += "volume " + str(L**3) + "\n"
s += "temperature 1.0\n"
s += "cutoff 1.0\n"
s += "boundary halo " + str(halo) + "\n\n"

s += "timestep " + str(dt) + "\n"
s += "steps " + str(N) + "\n"
s += "equilibration steps 0\n"
s += "scale temperature every 10\n"
s += "trajectory 0 100\n"
s += "stats every 100\n"
s += "stack size 100\n"
s += "print every " + str(thermo) + "\n\n"

s += "job time 1000000.0\n"
s += "close time 1.0\n\n"

s += "ensemble nvt dpdvv\n\n"
s += "finish\n"

print("Box size: %.1f | Timestep: %.3f | Num steps: %i" % (L, dt, N))
open("CONTROL", "w").write(s)
print("CONTROL file saved.")

