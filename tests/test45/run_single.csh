#!/bin/csh -f
#trace on
#---------------------------------------------------------------
#
#  Author A.Ivanchenko 27 May 2005
#  test45
#
#----------------------------------------------------------------

echo 'Start execute for target ' $TARGET 

mkdir -p  $TARGET
cd  $TARGET
rm -f r.out 

$G4MY/test45 $TEST45/$TARGET/run.mac >& r.out

cd ../
 
source $TEST45/plot.csh 
 




