#!/bin/bash

usage() { echo "$0:
  Script to extract eigenvalues from Gaussian output files
  to compute density of states
  Last change: 13/04/15
  Cmd ln args:
  1  base_directory, e.g. ~/Platinum/Plain
  2  cluster type, e.g. 9_10_9"
}

if [ $1 == "-h" ]; then
  usage
  exit 0
fi

base_dir=$1
cluster=$2

#
# create lists of eigenvalues for each spin
#
for i in {0..10}; do
  dir=$base_dir/Pt$cluster/S_$i
  >$dir/evals.out
  >$dir/evals_occ.out
  a=`cat $dir/Pt.out | grep "eigenvalues" | awk '{for(i=5; i<=NF; i++) print $i}'`
  b=`cat $dir/Pt.out | grep "occ" | awk '{for(i=5; i<=NF; i++) print $i}'`
  
  if [ -n "$a" ]; then      # if string length is non-zero
    echo "$a" >>$dir/evals.out
    echo "$b" >>$dir/evals_occ.out
    # if S=0, states are degenerate and must double evals
    if [ $i -eq "0" ]; then          
      echo "$a" >>$dir/evals.out
      echo "$b" >>$dir/evals_occ.out
    fi
    num_a=`cat $dir/evals.out | wc -l`
    num_b=`cat $dir/evals_occ.out | wc -l`
    echo "$num_a: e-vals for Pt$cluster, S=$i."
    echo "$num_b: occ e-vals for Pt$cluster, S=$i."
  else
    echo "No eigenvalues in the output file."
  fi
done

#
# paste eigenvalues for each spin in a table
#
if [ ! -d "$base_dir/Outfiles" ]; then
  mkdir $base_dir/Outfiles
fi
outdir=$base_dir/Outfiles/Evals
if [ ! -d $outdir ]; then
  mkdir $outdir
fi
file_a=evals_Pt${cluster}.out
file_b=evals_Pt${cluster}_occ.out
>$outdir/$file_a
>$outdir/$file_b

for i in {0..10}; do
  dir=$base_dir/Pt$cluster/S_$i
  paste -d " " $outdir/$file_a $dir/evals.out >$outdir/temp_a
  paste -d " " $outdir/$file_b $dir/evals_occ.out >$outdir/temp_b
  mv $outdir/temp_a $outdir/$file_a
  mv $outdir/temp_b $outdir/$file_b
  rm $dir/evals.out $dir/evals_occ.out
done

echo "Eigenvalues for Pt$cluster saved in $outdir/$file_a and $outdir/$file_b."

