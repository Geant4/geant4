#!/bin/sh -f
#
#  This script takes a gmake.log file
# produced by a reconstruction. It filters
# this file by removing lines found in 
# file in a specific .filter file.
#
# Usage :
#     UNIX> cd <g4install>/tests/tools/bin
#     UNIX> chmod u+x filter.sh
#     UNIX> ./filter.sh OSF1
#
# Some checks :
if [ -z "${G4INSTALL}" ] ; then
  echo "G4INSTALL not set. Execute setup file first !"
  exit
fi
if [ -z "${G4WORKDIR}" ] ; then
  echo "G4WORKDIR not set. Execute setup file first !"
  exit
fi
grep -v -F -f $G4INSTALL/tests/tools/bin/$1.filter $G4WORKDIR/stt/$1/gmake.log
