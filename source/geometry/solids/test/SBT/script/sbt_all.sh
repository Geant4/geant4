#!/bin/sh
# 
# This script runs all tests in SBT/geant4 directory 
# and stores outputs from each test in log/SBT.solid.out files.

#set -x

CURDIR=`pwd`

for SOLID_MACRO in `ls geant4/*.geant4`
do
  SOLID=`echo $SOLID_MACRO | sed 'skgeant4/kk' | sed 'sk.geant4kk'`
  echo "... Running SBT test for $SOLID" 
  $CURDIR/script/sbt.sh $SOLID
done  
  
  
