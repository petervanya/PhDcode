#!/bin/bash

nb=$1

if [ -z "$nb" ]; then
    echo -e "Usage:\n    parse_surf_ten.sh <nb>"
    exit 0
fi

echo "Extracting surface tensions for bin size: $nb"

dirn=Gamma_nb$nb
outname=gamma_raw.out
if [ ! -d "$dirn" ]; then mkdir $dirn; fi

for conf in `cat configs.out `; do 
    if [ ! -f "$conf/surf_nb$nb.out" ]; then
        echo "No surf_nb$nb.out found in $conf."
        exit 0
    fi
    gamma=`tail $conf/surf_nb$nb.out | tail -n 1 | awk '{print $3, $6}'`
    echo $conf $gamma
done >$dirn/$outname

echo "File $outname saved in $dirn/."
