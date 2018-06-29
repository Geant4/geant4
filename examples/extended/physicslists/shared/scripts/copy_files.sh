#!/bin/sh

# Script to copy classes from shared directory into include and src
# as it is required by GNUmake build.  Create the subdirs if they
# don't exist.
# Usage: copy_files.sh directoryName
#
# By I. Hrivnacova, IPN Orsay

DIRNAME=$1

# make subdirs if necessary
if [ ! -d include ]; then mkdir include ; fi
if [ ! -d src     ]; then mkdir src     ; fi
# copy files from $DIRNAME
cp -rp $DIRNAME/include/* include
cp -rp $DIRNAME/src/* src

echo "... copy_files.sh from $1 finished"

