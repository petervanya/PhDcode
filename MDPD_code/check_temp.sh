#!/bin/bash

for dirn in `cat configs.out`; do
    printf "%20s  " $dirn
    if [ -f "$dirn/OUTPUT" ]; then
        temp=`cat $dirn/OUTPUT | grep "temperature" -A 2 | \
            tail -n 1 | awk '{print $NF}'`
        echo "$temp"
    else
        echo "not"
    fi
done
cd ..

