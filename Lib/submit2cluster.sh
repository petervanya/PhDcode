#!/bin/bash

Usage(){
  echo "Usage: 
    submit2cluster.sh <prog> <file> <nc>
  
Script to submit Gaussian adsorption runs to the Cottrell

Arguments:
    <prog>        Programme to use, 'gaussian' or 'lammps'
    <file>        Absolute path to the file with NO extension

Options:
    <nc>          Number of cores"
exit 0
}


if [ "$1" = "-h" ]; then
    Usage
fi

# Comment due to DLMS not requiring input file arguments
if [ -z "$2" ]; then
    Usage
else
    filepath=$2
fi

if [ -z "$3" ]; then
    Ncores=16
else
    Ncores=$3
fi

# ======================= LEAVE THIS HERE!
#$ -cwd                                                                         
#$ -j y                                                                         
#$ -S /bin/bash                                                                 

# Gaussian commands                                                             
g09root="/home/Gaussian"
#  GAUSS_SCRDIR="/raid/home/Gaussian_scratch"
GAUSS_EXEDIR="/home/Gaussian/g09/bsd:/home/Gaussian/g09/private:/home/Gaussian/g09"
export g09root GAUSS_SCRDIR GAUSS_EXEDIR
. $g09root/g09/bsd/g09.profile
# =======================

filedir=$(dirname "$filepath")
filename=$(basename "$filepath")
cd $filedir

if [ "$1" == "gaussian" ]; then
    echo "Running Gaussian."
    progpath="/home/Gaussian/g09/g09"
    $progpath <$filepath.gjf >$filepath.out
elif [ "$1" == "lammps" ]; then
    echo "Running LAMMPS."
    progpath="/home/pv278/sw/bin/lmp_mpi"
    mpi_run="/opt/openmpi/bin/mpirun"
    time $mpi_run -np $Ncores $progpath <$filepath
elif [ "$1" == "dlms" ]; then
    echo "Running dl_meso_mpi."
    progpath="/home/pv278/sw/bin/dl_meso_mpi"
    mpi_run="/opt/openmpi/bin/mpirun"
    time $mpi_run -np $Ncores $progpath
elif [ "$1" == "dlms_serial" ]; then
    echo "Running dl_meso_serial."
    progpath="/home/pv278/sw/bin/dl_meso_serial"
    time $progpath
else
    Usage
fi


