#!/bin/sh

# Script to copy classes from shared directory into include and src
# as it is required by GNUmake build.
# Usage: copy_files.sh directoryName
#
# By I. Hrivnacova, IPN Orsay

DIRNAME=$1

cp -rp $DIRNAME/include include
cp -rp $DIRNAME/src src

echo "... copy_files.sh from $1 finished"

