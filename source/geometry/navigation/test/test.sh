#!/bin/sh
echo "Run on `hostname`, which is a `uname -a` machine" | tee -a test.out
echo `date` >> test.out
G4BIN="$G4WORKDIR/bin/$G4SYSTEM"

for i in *.cc
do
  j=`basename $i .cc`
  make G4TARGET=$j
  echo Test output for $j... >> test.out
  $G4BIN/$j >> test.out 2>&1;
done


