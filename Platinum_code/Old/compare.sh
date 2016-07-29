#!/bin/bash

# USAGE
# =====
# Script to get the binding energy of water and carbon
# 1st cmd line arg -- number of Pt atoms
# =====

cd ../Oxygen/Pt$1/
a=`cat Pt$1.out | grep "SCF Done" | tail -1 | awk '{print $5}'`
b=`cat Pt$1_plain.out | grep "SCF Done" | awk '{print $5}'`
#c=-76.4087650029   # water
c=-75.0606231181   # oxygen

echo "X+oxygen = $a"
echo "X = $b"
echo "oxygen = $c"

echo "$a-($b $c)" | bc -l | awk '{print "Energy difference = " $1*27.2 " eV"}'

