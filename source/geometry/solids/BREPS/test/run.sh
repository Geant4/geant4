#!/usr/local/bin/bash
echo "Running on `hostname`, which is a `uname -a` machine" 
k=$TESTTARGET
# if ( ! -f $G4SYSTEM ) ln -fs ../../../../bin/$G4SYSTEM .
for i in CurveTest.cc G4*.cc
do
  j=`basename $i .cc`
  TESTTARGET=$j
  export TESTTARGET
  gmake
  echo Test output for $TESTTARGET... 
  #$G4SYSTEM/$TESTTARGET 
  #../../../../bin/$G4SYSTEM/$TESTTARGET 
  ./bin/$G4SYSTEM/$TESTTARGET 
done

# Now test STEPtest
export TESTTARGET=STEPTest
gmake
echo Test outputs for $TESTTARGET... 
for j in 1 2 3 4 5 6 7 8 9 
do
  echo $j | ./bin/$G4SYSTEM/$TESTTARGET
done

TESTTARGET=$k
export TESTTARGET
