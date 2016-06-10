#!/bin/bash
#
# Usage:
# check_tabs.sh 

#set -x

CURDIR=`pwd`

for SOURCE in `ls *.cc include/*.hh src/*.cc`
do
  RESULT=`cat $SOURCE | sed 's/\t/TAB_TO_BE_REMOVED/g' | grep TAB_TO_BE_REMOVED` 
  if [ ! "$RESULT" = "" ]; 
  then 
    echo "TAB detected in $SOURCE"
    #cat $SOURCE | sed 's/\t/TAB_TO_BE_REMOVED/g' > $CURDIR/"$DIR"_tabs/$SOURCE
  fi
done 

cd $CURDIR
