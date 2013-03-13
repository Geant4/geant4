#!/bin/bash
#
# Usage:
# replace_tabs.sh 
#
# Script for fixing Geant4 examples coding guideline 3.1 (tabulations).
# By I. Hrivnacova, IPN Orsay

#set -x

CURDIR=`pwd`

for SOURCE in `ls *.cc include/*.hh src/*.cc`
do
  RESULT=`cat $SOURCE | sed 's/\t/TAB_TO_BE_REMOVED/g' | grep TAB_TO_BE_REMOVED` 
  if [ ! "$RESULT" = "" ]; 
  then 
    echo "Updating file $SOURCE ... "
    cat $SOURCE | sed --posix -e 's/\t/        /g' > $SOURCE.1
    mv $SOURCE.1 $SOURCE
  fi  
done 

cd $CURDIR
