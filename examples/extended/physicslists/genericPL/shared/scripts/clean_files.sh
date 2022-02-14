#!/bin/sh

# Script to remove classes copied from shared directory into include and src
# via copy_files.sh
#
# By I. Hrivnacova, IPN Orsay

DIRNAME=$1

for FILE in `ls $DIRNAME/include`; do
  rm -f include/$FILE
done
for FILE in `ls $DIRNAME/src`; do
  rm -f src/$FILE
done

echo "... clean_files.sh from $1 finished"

