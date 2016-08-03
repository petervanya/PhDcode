#!/bin/bash
# Temporary production of input.yaml file
# 28/07/16

fname="input.yaml"

cat >$fname <<EOF
# Methods:
# 1. aii = 25 kT, da ~chi from FH
# 2. aii based on compressibility, da = sqrt(a11 a22) + ~chi
# 3. from Travis, JCP, 2007

dpd-density: 3.0
box-size: 10.0            # in nm
temperature: 300.0        # in K
dt: 0.05
run-time: 10              # in ns
water-meoh-ratio: 1.0     # number ratio of water to meoh molecules
method: 1
waters-in-bead: 2
meohs-in-bead: 1
EOF
echo "Parameter file saved in $fname"
