#!/bin/sh
# 
# This script must be run from the top level directory where is SBT installed
##############################################################################

if [ $# -ne 2 ]
then
  echo
  echo "Usage: $0 sbtexecutable solidname"
  echo
  exit 1
fi

echo `pwd`

export EXE=$1
shift
export SOLID=$1
shift

export INPUT=geant4/${SOLID}.geant4
export LOG=log/SBT.${SOLID}.log

export COUNTERR="awk -f script/counterr.awk log/${SOLID}.*.log"
export COUNTVOXELERR="awk -f script/countvoxelerr.awk log/${SOLID}v.*.log"

rm -f ${LOG}
rm -f sbt.out

time ${EXE} < ${INPUT} 2>&1 > ${LOG} 2>&1 > sbt.out
time ${COUNTERR} 2>&1 >> sbt.out
time ${COUNTVOXELERR} 2>&1 >> sbt.out
