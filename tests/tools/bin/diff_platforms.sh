#!/bin/sh
# $Id: diff_platforms.sh,v 1.1 1999-01-08 16:36:05 gunter Exp $
# Produces differences of *.out files between platforms.
# The result is to standard output.
# Usage:         diff_platforms.sh <platform1> <platform2>
# E.g.:          diff_platforms.sh DEC6-AFS SUN-AFS
# Suggested use: diff_platforms.sh DEC6-AFS SUN-AFS > diff.lis

if [ $# -ne 2 ]
then
  echo 'Requires 2 arguments, <platform1> <platform2>.'
  exit
fi

basedir=$G4WORKDIR/..
dir1=$basedir/$1/stt/$1
dir2=$basedir/$2/stt/$2

for i in $dir1/*.out
do
  j=`basename $i`
  echo \\n\\nDiffing $j for $1 and $2...\\n\\n
  diff $dir1/$j $dir2/$j
done

exit
