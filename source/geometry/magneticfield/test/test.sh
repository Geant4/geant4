#!/bin/sh
#
# A simple script to run all the tests in this directory and check
# their results against the expected (previous) results
#
# $Id$
# $Name: geant4-09-02-ref-05 $
#

echo "Running on `hostname`, which is a `uname -a` machine" 
host=`hostname`

#
MAKE=make
BINDIR=$G4BIN/$G4SYSTEM/
#xport BINDIR=$G4BIN_GMAKE/$G4SYSTEM/
#INDIR=$G4WORKDIR/bin/$G4SYSTEM/

target=testPropagateMagField
echo  "Compiling $target ... "
$MAKE G4TARGET=$target   || exit
if [[ ! -x $BINDIR/$target ]] ; then
  echo "Could not find executable $target in directory $BINDIR" 
  exit 1
fi
echo  "Executing $target ..."
for n in 1 2 3 4 5 6 7 8 9 10 11 12 13
do
  echo "Executing with stepper choice $n ..."
  $BINDIR/$target $n  > $target.newout$n \
		     2> $target.newerr$n
  echo  ".. difference from expected output: "
  diff -wb $target.out$n $target.newout$n
  sleep 1;
  echo  ".. difference from expected error: "
  diff -wb $target.err$n $target.newerr$n
  sleep 1;
  echo  " "
done

target=testProElectroMagField
echo  "Compiling $target ... "
$MAKE G4TARGET=$target  || exit
echo  "Executing $target ..."
for n in 1 2 3 4 8
do
  echo "Executing with stepper choice $n ..."
  $BINDIR/$target $n  > $target.newout$n \
		     2> $target.newerr$n
  echo  ".. difference from expected output: "
  diff -wb $target.out$n $target.newout$n
  sleep 1;
  echo  ".. difference from expected error: "
  diff -wb $target.err$n $target.newerr$n
  sleep 1;
  echo  " "
done

for i in *Spin.cc
do
  target=`basename $i .cc`
  echo  "Compiling $target ... "
  $MAKE G4TARGET=$target   || exit
  echo  "Executing $target ..."
  for n in  4 3 2 1 0
  do
    echo "Executing with stepper choice $n ..."
    $BINDIR/$target $n  > $target.newout$n \
		       2> $target.newerr$n
    echo  ".. difference from expected output: "
    diff -wb $target.out$n $target.newout$n
    sleep 1;
    echo  ".. difference from expected error: "
    diff -wb $target.err$n $target.newerr$n
    sleep 1;
    echo  " "
  done
done

exit

for i in *Spin.cc
do
  target=`basename $i .cc`
  echo  "Compiling $target ... "
  gmake G4TARGET=$target || exit
  echo  "Executing $target ..."
  $BINDIR/$target $n  > $target.newout$n \
		     2> $target.newerr$n
  echo  ".. difference from expected output: "
  diff -wb $target.out $target.newout
  sleep 1;
  echo  ".. difference from expected error: "
  diff -wb $target.err $target.newerr
  sleep 1;
  echo  " "
done

exit
