#!/bin/sh -f
#
#  This script filters a :
#    $G4WORKDIR/stt/$G4SYSTEM/gmake.log 
# file resulting from a reconstruction.
# It does a specific awk over the filtered output and
# produces a : 
#     $G4WORKDIR/stt/$G4SYSTEM/problems
# output file.
#  A <g4system>.awk file should be provided
# for each <g4system>. See OSF1.awk, HP-aCC.awk
# for examples.
#
#  Usage :
#     UNIX> cd <g4install>/tests/tools/bin
#     UNIX> chmod u+x analyse.sh
#     UNIX> ./analyse.sh $G4SYSTEM
#    (UNIX> ./analyse.sh OSF1)
#
file=$G4WORKDIR/stt/$1/problems
$G4INSTALL/tests/tools/bin/filter.sh $1 | awk -f $G4INSTALL/tests/tools/bin/$1.awk > $file
grep Error_in $file
grep Warning_in $file
grep Compiler_crash_in $file
if [ -s $file ] ; then
  echo "See file $file for details."
fi
#