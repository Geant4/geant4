#!/bin/sh
#
# A simple script to run all the tests in this directory and check
# their results against the expected (previous) results
#
# $Id: test.sh,v 1.2 2000-11-21 13:00:04 gcosmo Exp $
# $Name: not supported by cvs2svn $
#

echo "Running on `hostname`, which is a `uname -a` machine" 
host=`hostname`

target=testPropagateMagField
echo  "Compiling $target ... "
gmake G4TARGET=$target
echo  "Executing $target ..."
for n in 1 2 3 4 5 6 7 8
do
  echo "Executing with stepper choice $n .. \c"
  $G4WORKDIR/bin/$G4SYSTEM/$target $n > $target.newout$n \
                                     2> $target.newerr$n
  echo  ".. difference from expected output: "
  diff -wb $target.out$n $target.newout$n
  sleep 1;
  echo  ".. difference from expected error: "
  diff -wb $target.err$n $target.newerr$n
  sleep 1;
  echo  " "
done

#exit

for i in *Spin.cc testProElectroMagField.cc
do
  target=`basename $i .cc`
  echo  "Compiling $target ... "
  gmake G4TARGET=$target
  echo  "Executing $target ..."
  $G4WORKDIR/bin/$G4SYSTEM/$target > $target.newout \
                                  2> $target.newerr
  echo  ".. difference from expected output: "
  diff -wb $target.out $target.newout
  sleep 1;
  echo  ".. difference from expected error: "
  diff -wb $target.err $target.newerr
  sleep 1;
  echo  " "
done

exit
