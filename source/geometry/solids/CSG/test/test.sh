#!/usr/local/bin/bash
echo `date` >> test.out
k=$TESTTARGET
for i in *.cc
do
  j=`basename $i .cc`
  TESTTARGET=$j
  export TESTTARGET
  gmake
  echo Test output for $TESTTARGET... >>test.out
  $G4SYSTEM/$TESTTARGET >>test.out 2>&1;
done
TESTTARGET=$k
export TESTTARGET


