#!/bin/ksh

if [ $# -lt 3 ]
then
   echo Usage: streamedit.sh oldstring newstring files
   exit 1
fi

echo Replacing all occurrences of $1 with $2 \
     in the following files:

a=$1; shift;
b=$1; shift;
for i
do
  echo $i
  cp $i $i.bak
  perl -pe "s#$a#$b#g;" <$i >$i.new
  cp $i.new $i
done

