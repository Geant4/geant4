#!/bin/sh
# $Id: diff_ref_ref+_outputs.sh,v 1.2 1999-07-20 09:55:21 stesting Exp $
# Produces difference of prod and dev .out files.
# The result is in both prod's and dev's  $G4WORKDIR/stt/$G4SYSTEM/ in
#   a file called diff_prod_dev_outputs.lis.

prodworkdir=`echo $G4WORKDIR | sed s/dev/prod/`
devworkdir=`echo $prodworkdir | sed s/prod/dev/`
proddir=$prodworkdir/stt/$G4SYSTEM
devdir=$devworkdir/stt/$G4SYSTEM

lisprodfile=$proddir/diff_prod_dev_outputs.lis
lisdevfile=$devdir/diff_prod_dev_outputs.lis
rm -f $lisprodfile
rm -f $lisdevfile
touch $lisprodfile
touch $lisdevfile

for i in $proddir/*.out
do
  j=`basename $i`
  diff $proddir/$j $devdir/$j >>$lisprodfile 2>&1
done
cp $lisprodfile $lisdevfile

echo Differences are in \$G4WORKDIR/stt/$G4SYSTEM/diff_prod_dev_outputs.lis.
