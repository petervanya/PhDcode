#!/bin/bash
# =====
# Summarised table of each spin config for a Pt cluster from Gaussian outfiles
# Spin | converged? | Energy
# Cmd ln args
# 1 -- molecule, e.g. Water
# 2 -- cluster, e.g. 12_7
# 3 -- position of water on Pt cluster 1-3
# =====

dir=~/Adsorption/$1
cluster=Pt$2
if [ $1 == "Water" ]; then
	cluster=Pt$2/Eta_$3
fi
cd $dir/$cluster

#echo -e "# Spin \t Conv \t Num. cycles \t Eerr \t Runtime"

for i in {0..10}
do
  if [ `tail -n 1 S_$i/Pt.out | awk '{print $1}'` == "Normal" ]
	then 
		conv="Yes"
    #ZPE=`cat S_$i/Pt.out | grep "Sum of electronic and zero-point" | awk '{print $NF}'`
    ZPE=`cat S_$i/Pt.out | grep "SCF Done" | tail -n 1 | awk '{print $5}'`
  else
		conv="No"
	  ZPE="NA"
  fi

	echo -e $i "\t" $conv "\t" $ZPE
done  
