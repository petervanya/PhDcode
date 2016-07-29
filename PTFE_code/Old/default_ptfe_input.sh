#!/bin/bash
# Default input.yaml serving for creation of LAMMPS data file to simulate Nafion w carbon electrodes
# 08/06/16

cat >input.yaml <<EOF
# ===== Bead types:
# * A, B, C
# * W: water bead 6 H2O
# * E: electrodes from carbon/quartz/silica
# * P: platinum
# =====
box-size:          40           # DPD units, 1 = 8.14 AA
temperature:       1.0          # Units of kB T
topology:          (A 3, [B 1], [C 1])15
mono-per-chain:    15
water-uptake:      9            # water uptake, number of H2O/SO3H
gamma:             4.5          # DPD drag coefficient

chi-params:                     # DPD Flory-Huggins params (Wu, EES, 2008)
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

bond-coeffs:                    # spring const param k_ij
    A A: 4.0
    A B: 4.0
    A C: 4.0

equilibrium-dist:  0.1          # r0 in k(r-r0)^2
# spring-const:     4.0

electrodes:
    width:         5.0          # width of one electrode, same on the other side
    Pt-amount:     0.0          # number of Pt beads per number of carbon beads
EOF
echo "Parameter file saved in input.yaml"
