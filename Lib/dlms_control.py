#!/usr/bin/env python
"""Usage:
    dlms_control.py [--L <L> --dt <dt> --steps <n> --startstep <st>]
                    [--thermo <th> --halo <h> --eq <Neq>]

Generate DL_MESO control file.

Options:
    --L <L>            Box size, either 'L' or 'Lx Ly Lz' [default: 40.0]
    --dt <dt>          Timestep [default: 0.05]
    --steps <n>        Number of simulation steps [default: 10000]
    --startstep <st>   Start timestep for dumping frames [default: 0]
    --eq <Neq>         Equilibration steps [default: 0]
    --thermo <th>      Print every [default: 100]
    --halo <h>         Boundary halo, like neighbor [default: 2.5]

pv278@cam.ac.uk, 06/06/16
"""
from docopt import docopt
import sys


args = docopt(__doc__)
sL = args["--L"].split()
if len(sL) == 1:
    L = list(map(float, sL)) * 3
elif len(sL) == 3:
    L = list(map(float, sL))
else:
    sys.exit("<L> should have size 1 or 3.")
dt = float(args["--dt"])
N = int(args["--steps"])
startstep = int(args["--startstep"])
thermo = int(args["--thermo"])
halo = float(args["--halo"])
eqsteps = int(args["--eq"])


s = "pokus\n\n"

s += "volume %.1f %.1f %.1f\n" % (L[0], L[1], L[2])
s += "temperature 1.0\n"
s += "cutoff 1.0\n"
s += "boundary halo %.1f\n\n" % halo

s += "timestep %.3f\n" % dt
s += "steps %i\n" % N
s += "equilibration steps %i\n" % eqsteps
s += "scale temperature every 10\n"
s += "trajectory %i 100\n" % startstep
s += "stats every 100\n"
s += "stack size 100\n"
s += "print every %i\n\n" % thermo

s += "job time 1000000.0\n"
s += "close time 1.0\n\n"

s += "ensemble nvt dpdvv\n\n"
s += "finish\n"

print("Box: [ %.1f  %.1f  %.1f] | Timestep: %.3f | Num. steps: %i" % \
        (L[0], L[1], L[2], dt, N))
open("CONTROL", "w").write(s)
print("CONTROL file saved.")


