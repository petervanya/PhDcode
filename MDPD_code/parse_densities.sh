#!/bin/bash

dirn=Density
outname=density_raw.out
if [ ! -d "$dirn" ]; then mkdir $dirn; fi

for conf in `cat configs.out `; do 
    rho=`cat $conf/density.out | grep "density:" | awk '{print $NF}'`
    echo $conf $rho
done >$dirn/$outname
echo "File $outname saved in $dirn/."
