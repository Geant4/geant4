#!/bin/bash
#
# Usage:
# check_separator.sh 

#set -x

CURDIR=`pwd`

#The correct separator with 80 characters
#SEPARATOR="//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......"

#Still tolerated separator with 78 characters
SEPARATOR="//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...."

for SOURCE in `ls src/*.cc`
do
  LINES=`cat $SOURCE | grep $SEPARATOR | wc -l `  
  RESULT=`echo "$LINES > 0" | bc`
  if [ $RESULT -ne 1 ]; 
  then 
    echo "NO SEPARATOR found in $SOURCE"
  fi
done 

cd $CURDIR
