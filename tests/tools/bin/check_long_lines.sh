#!/bin/bash
#
# Usage:
# check_long_lines.sh [maxColumn]
#
# Script for checking Geant4 examples coding guideline 3.1.
# By I. Hrivnacova, IPN Orsay

#set -x

CURDIR=`pwd`

MAX=90
if [ $# -eq 1 ]; then
  MAX=$1
fi  

for SOURCE in `ls *.cc include/*.hh src/*.cc`
do
  LINE_LENGTH=`wc -L $SOURCE | sed sY" $SOURCE"YYg` 
  RESULT=`echo "$LINE_LENGTH > $MAX" | bc`
  if [ $RESULT -ne 0 ]; 
  then 
    echo "LONG LINE ($LINE_LENGTH characters) detected in $SOURCE"
  fi
done 

cd $CURDIR
