#!/usr/local/bin/bash
echo "Run on `hostname`, which is a `uname -a` machine" | tee -a test.out
echo `date` >> test.out
k=$TESTTARGET
ln -fs ../../../../bin/$G4SYSTEM .
for i in *.cc
do
  j=`basename $i .cc`
  TESTTARGET=$j
  export TESTTARGET
  gmake
  echo Test output for $TESTTARGET... >>test.out
  $G4SYSTEM/$TESTTARGET >>test.out 2>&1;
  # ../../../../bin/$G4SYSTEM/$TESTTARGET >>test.out 2>&1;
done
TESTTARGET=$k
export TESTTARGET


