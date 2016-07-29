#!/bin/bash

# Summarised table of of water on Pt adsorption runs
# for each position (top, bridge, fcc) and each spin


dir=~/Platinum/Water
cluster=Pt$1
if [ ! -d $dir/$cluster ]; then
  echo "$cluster does not exist"
  exit 1
fi
cd $dir/$cluster

echo -e "Eta \t Spin \t Succ? \t Reason \t Steps \t MaxSteps \t E \t" \
         "Runtime d:h:m:s \t Date"

for j in {1..3}; do
  for i in {0..10}; do
    filecheck="Eta_$j/S_$i/Pt.out"
    if [ ! -e $filecheck ]; then
      echo "$filecheck: file does not exist."
      continue
    fi
    if [ `cat $filecheck | wc -l` == "0" ]; then
      echo "$filecheck: empty file."
      continue
    fi

    cd Eta_$j/S_$i
    spin=`grep "^ Charge" Pt.out | head -n 1 | awk '{print substr($NF,length($NF)-1,length($NF))}'`
    spin=`echo "($spin-1)/2" | bc`
    if [ `tail -n 1 Pt.out | awk '{print $1}'` == "Normal" ]; then 
      succ="Yes"
      reason="NA\t"
      E=`cat Pt.out | grep "cycles$" | tail -1 | awk '{print $5}'`
    else 
      succ="No"
      E="NA\t"
      str=`tail -n 4 Pt.out | head -n 1 | awk '{print $1}'`
      if [ $str == "Convergence" ]; then
        reason="Conv failure"
      elif [ $str == "Error" ]; then
        reason="Ran out of steps"
      else 
        reason="NA"
      fi

    fi

    steps=`cat Pt.out | grep "^ Step number" | tail -1 | awk '{print $3}'`
    maxsteps=`cat Pt.out | grep "^ Step number" | tail -1 | awk '{print $NF}'`
    runtime=`cat Pt.out | grep "^ Job cpu" | awk '{print $4":"$6":"$8":"$10}'`
    datetime=`cat Pt.out | grep " termination " | tail -1 | awk '{print $(NF-3),$(NF-2),$NF,$(NF-1)}'`
    
    echo -e $j "\t" $spin "\t" $succ "\t" $reason "\t" $steps "\t" $maxsteps "\t" $E "\t" $runtime "\t" $datetime
    cd ../..
  done
done

