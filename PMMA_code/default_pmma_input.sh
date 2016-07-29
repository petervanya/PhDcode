#!/bin/bash

file=input.yaml

cat >$file << EOF
# Input file for the simulation of PMMA with solvent (water or methanol)
# NOW OBSOLETE, MODIFY FOR DPD UNITS

solvent-vol:    0.2
n:              100              # PMMA polymerisation
temperature:    300
box-size:       40               # in DPD units
gamma:          4.5              # as in Groot_JCP_1997

ksi-params:
    1 2:        17.1             # PMMA-water

bond-coeffs:
    1 2:        4.0
EOF
echo "Input file saved in $file"
