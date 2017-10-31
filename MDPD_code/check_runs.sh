#!/bin/bash

for dirn in `cat configs.out`; do
   printf "%20s  " $dirn 
   if [ -f "$dirn/OUTPUT" ]; then
       line=`tail -n 3 $dirn/OUTPUT | head -n 1 | awk '{print $1, " ", $NF}'`
       printf "$line\n"
   else
       echo "not"
   fi
done

