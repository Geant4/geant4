#!/usr/local/bin/bash
#
# A simple script to run all the tests in this directory and check
#  their results against the expected (previous) results
#
# $ Id: $
# $ Name : $
#
# Created:  12 May 99 J. Apostolakis: starting from P.Kent's test.sh
#                                     changed output to STDout only        
# Modified: 21 May 99 J. Apostolakis: check the results

echo "Running on `hostname`, which is a `uname -a` machine" 
host=`hostname`
k=$TESTTARGET
# if ( ! -f bin ) ln -fs ../../../../bin . 
for i in CurveTest.cc G4*.cc
do
  target=`basename $i .cc`
  TESTTARGET=$target
  export TESTTARGET
  echo  "Compiling $target ... "
  gmake
  echo  "Executing $target ... "
  ./bin/$G4SYSTEM/$target > $target.newout-$host
  echo  "Difference from expected output: "
  diff $target.out $target.newout-$host
  echo  " "
done

exit

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
