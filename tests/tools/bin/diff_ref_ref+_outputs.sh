#!/bin/sh
# $Id: diff_ref_ref+_outputs.sh,v 1.1 1999-01-08 16:36:05 gunter Exp $
# Produces difference of ref and ref+ .out files.
# The result is in both ref's and ref+'s  $G4WORKDIR/stt/$G4SYSTEM/ in
#   a file called diff_ref_ref+_outputs.lis.

refworkdir=`echo $G4WORKDIR | sed s/ref+/ref/`
refplusworkdir=`echo $refworkdir | sed s/ref/ref+/`
refdir=$refworkdir/stt/$G4SYSTEM
refplusdir=$refplusworkdir/stt/$G4SYSTEM

lisfile=$refdir/diff_ref_ref+_outputs.lis
lisplusfile=$refplusdir/diff_ref_ref+_outputs.lis
rm -f $lisfile
rm -f $lisplusfile
touch $lisfile
touch $lisplusfile

for i in $refdir/*.out
do
  j=`basename $i`
  diff $refdir/$j $refplusdir/$j >>$lisfile 2>&1
done
cp $lisfile $lisplusfile

echo Differences are in \$G4WORKDIR/stt/$G4SYSTEM/diff_ref_ref+_outputs.lis.
