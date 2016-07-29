#!/bin/bash

# USAGE
# =====
# -- 1st cmd line arg: config, e.g. 5_10_5
# -- 2nd cmd line arg: spin
# =====

# Generic Cottrell header
# =======================
#$ -cwd                                                                         
#$ -j y                                                                         
#$ -S /bin/bash                                                                 

# Gaussian commands                                                             
g09root="/home/Gaussian"
GAUSS_SCRDIR="/state/partition1/Gaussian_scratch"
GAUSS_EXEDIR="/home/Gaussian/g09/bsd:/home/Gaussian/g09/private:/home/Gaussian/g09"
export g09root GAUSS_SCRDIR GAUSS_EXEDIR
# =======================

config=$1          # e.g. 5_10_5
spin=$2            # e.g. 1

cd /home/pv278/Adsorption/Plain/LANL2MB/Pt$config/S_$spin
/home/Gaussian/g09/g09 < Pt.gjf > Pt.out
