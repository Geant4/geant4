#!/bin/sh -f
#
#  To help finding source file.
#
#  Usage :
#     UNIX> cd <g4install>/tests/tools/bin
#     UNIX> chmod u+x search.sh
#     UNIX> ./search.sh G4ElScatterer.cc
#
find $G4INSTALL/source -name "*$1*" -print
