#!/bin/sh
#
# $Id: sbt.sh,v 1.2 2008-01-25 17:09:50 ivana Exp $
# 
# This script must be run from the top level directory where is SBT installed;
# it calls the SBT program with the selected geant4 macro in SBT/geant4 directory.
# The argument has to specify the solid name.

if [ $# -ne 1 ]
then
  echo
  echo "Usage: sbt.sh solidname"
  echo
  exit 1
fi

SOLID=$1
INPUT=geant4/${SOLID}.geant4
OUTPUT=log/SBT.${SOLID}.out
LOG=log/SBT.${SOLID}.log

COUNTERR="awk -f script/counterr.awk log/${SOLID}.*.log"
COUNTVOXELERR="awk -f script/countvoxelerr.awk log/${SOLID}v.*.log"

rm -f ${LOG}
rm -f ${OUTPUT}
rm -f tmp.out

# Run the test
#
{ time SBT < ${INPUT} >& ${LOG}; } >& tmp.out

# Print the time and error statistics
#
echo "SBT test for $SOLID"  >> ${OUTPUT}
echo >> ${OUTPUT}

# Run tests
if [ "`ls log/${SOLID}.*.log 2> /dev/null`" != "" ]
then 
  ${COUNTERR} >> ${OUTPUT}
else
  echo "No run test output"  >> ${OUTPUT}   
fi

# Voxel tests
if [ "`ls log/${SOLID}v.*.log 2> /dev/null`" != "" ]
then 
  ${COUNTVOXELERR} >> ${OUTPUT}
else
  echo "No voxel test output" >> ${OUTPUT}  
fi

echo >> ${OUTPUT}
echo "Time:"  >> ${OUTPUT}
cat tmp.out >> ${OUTPUT}
rm -f tmp.out
