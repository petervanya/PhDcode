#!/bin/bash

Usage(){
  echo "Usage: 
    $(basename $0) <sup>
    $(basename $0) pipe <user>
      
Check occupancy of Cottrell servers
  
Arguments:
    jae   Check JAE servers
    pdb   Check PDB servers
    pipe  Find own jobs waiting in the queue
  
pv278@cam.ac.uk, 16/05/15"
  exit 0
}

if [ "$1" == "-h" ]; then
  Usage
fi

RED='\033[00;31m'
GREEN='\033[00;32m'
RESTORE='\033[0m'

sup=$1
user=$2

if [ "`which qstat`" == "" ]; then
  echo "No qstat on this computer."
  exit 0
fi


if [ -z $sup ]; then
  echo -e "Usage: \n    occupancy.sh <sup>"
  exit
elif [ $sup == "jae" ]; then
  for i in {0..9}; do
    line=`qstat -f -u '*' | grep $sup.q@compute-0-$i -A 1 | tail -1`
    num=`echo $line | awk '{print $1}'`
    if [[ $num =~ ^[0-9]+$ ]] ; then
       person=`echo $line | awk '{print $4}'`
       cores=`qstat -f -u '*' | grep $sup.q@compute-0-$i | awk '{print $3}' | cut -d'/' -f2-3`
       echo -e $i ": ${RED}occupied${RESTORE}, $person, $cores"
    else
       echo -e $i ": ${GREEN}free${RESTORE}"
    fi
  done
elif [ $sup == "pdb" ]; then
  for i in {10..19}; do
    line=`qstat -f -u '*' | grep $sup.q@compute-0-$i -A 1 | tail -1`
    num=`echo $line | awk '{print $1}'`
    if [[ $num =~ ^[0-9]+$ ]] ; then
       person=`echo $line | awk '{print $4}'`
       cores=`qstat -f -u '*' | grep $sup.q@compute-0-$i | awk '{print $3}' | cut -d'/' -f2-3`
       echo -e $i ": ${RED}occupied${RESTORE}, $person, $cores"
    else
       echo -e $i ": ${GREEN}free${RESTORE}"
    fi
  done
elif [ $sup == "new" ]; then
  for i in {20..27}; do
    line=`qstat -f -u '*' | grep $sup.q@compute-0-$i -A 1 | tail -1`
    num=`echo $line | awk '{print $1}'`
    if [[ $num =~ ^[0-9]+$ ]] ; then
       person=`echo $line | awk '{print $4}'`
       cores=`qstat -f -u '*' | grep $sup.q@compute-0-$i | awk '{print $3}' | cut -d'/' -f2-3`
       echo -e $i ": ${RED}occupied${RESTORE}, $person, $cores"
    else
       echo -e $i ": ${GREEN}free${RESTORE}"
    fi
  done
elif [ $sup == "pipe" ]; then
  qstat -f -u '*' | grep "qw" | grep $user
  echo "`qstat -f -u '*' | grep "qw" | grep $user | wc -l` jobs of $user waiting"
else
  echo "Wrong argument, choose jae or pdb."
  exit 
fi
