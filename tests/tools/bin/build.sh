#!/bin/sh -f
#
#  This script spawn the reconstruction of libraries or tests
# from $G4INSTALL/tests/tools/bin directory. The output of the job
# is put in $G4WORKDIR/stt/$G4SYSTEM directory and is gmake.log.
#  Aliases g4tail, g4filter, defined by the setup.csh(.sh) file
# permits to follow the job, and to filter output.
#
#  Usage :
#     UNIX> cd <g4install>/tests/tools/bin
#     UNIX> chmod u+x build.sh
#     UNIX> ./build.sh all         - builds all libraries and builds all
#                                    test executables.
#     UNIX> ./build.sh lib all     - builds all libraries.
#     UNIX> ./build.sh lib <category>  - builds a specific category library.
#     UNIX> ./build.sh test all    - builds all test executables.
#     UNIX> ./build.sh test <test> - builds a specific test exectable.
#  And with csh :
#      csh> source setup.csh
#      csh> g4tail
#    ( csh> g4filter)
#
#  When env variables are correctly setted, this script could be 
# executed from everywhere.
#
# Some checks :
if [ -z "${MAKECMD}" ] ; then
MAKECMD='gmake -k'
fi
# Some checks :
if [ -z "${G4INSTALL}" ] ; then
  echo "G4INSTALL not set. Execute setup file first !"
  exit
fi
if [ -z "${G4SYSTEM}" ] ; then
  echo "G4SYSTEM not set. Execute setup file first !"
  exit
fi
if [ -z "${G4WORKDIR}" ] ; then
  echo "G4WORKDIR not set. Execute setup file first !"
  exit
fi
#
# Make $G4WORKDIR/stt directory :
dir=$G4WORKDIR/stt
if [ -d $dir ] ; then
 echo "" > /dev/null
else 
 echo "$dir does not exist."
 echo "This is now an error.  Create a directory or a symbolic link"
 echo "  to a directory by hand or in a calling script."
fi
#
# Make $G4WORKDIR/stt/$G4SYSTEM directory :
dir=$G4WORKDIR/stt/$G4SYSTEM
if [ -d $dir ] ; then
 echo "" > /dev/null
else 
 mkdir $dir
 echo "$dir created."
fi
#
# How to check if some argument given ?
#if [ -z $1 ] ; then 
#echo "One argument expected : libs or tests"
#exit
#fi
#
if [ "$1" = "lib" ] ; then 
if [ "$2" = "all" ] ; then 
# Go in source and spawn libs reconstruction :
cd $G4INSTALL/source
(date;time $MAKECMD;date) > $dir/gmake.log 2>&1
exit
else
# Go in source/$2 and spawn category reconstruction :
cd $G4INSTALL/source/$2
(date;time $MAKECMD;date) > $dir/gmake.log 2>&1
exit
fi
else
#
if [ "$1" = "test201" ] ; then 
# Go in tests/test201 and spawn reconstruction :
##CPPVERBOSE=1;export CPPVERBOSE
cd $G4INSTALL/tests/test201
(date;time $MAKECMD;date) > $dir/gmake.log 2>&1
else
# Do not use visualization :
G4VIS_USE_DAWN=;export G4VIS_USE_DAWN
G4VIS_USE_DAWNFILE=;export G4VIS_USE_DAWNFILE
G4VIS_USE_OPACS=;export G4VIS_USE_OPACS
G4VIS_USE_OPENGLX=;export G4VIS_USE_OPENGLX
G4VIS_USE_OPENGLXM=;export G4VIS_USE_OPENGLXM
G4VIS_USE_OIX=;export G4VIS_USE_OIX
G4VIS_USE_OIWIN32=;export G4VIS_USE_OIWIN32
G4VIS_USE_RAYX=;export G4VIS_USE_RAYX
G4VIS_USE_VRML=;export G4VIS_USE_VRML
G4VIS_USE_VRMLFILE=;export G4VIS_USE_VRMLFILE
#
if [ "$1" = "test" ] ; then 
if [ "$2" = "all" ] ; then 
# Go in tests and spawn reconstruction :
cd $G4INSTALL/tests
##CPPVERBOSE=1;export CPPVERBOSE
(date;time $MAKECMD;date) > $dir/gmake.log 2>&1
exit
else
# Go in tests/$2 and spawn reconstruction :
##CPPVERBOSE=1;export CPPVERBOSE
cd $G4INSTALL/tests/$2
(date;$MAKECMD clean_bin;time $MAKECMD;date) > $dir/gmake.log 2>&1
fi
else
#
if [ "$1" = "all" ] ; then 
(date;cd $G4INSTALL/source;time $MAKECMD;\
cd $G4INSTALL/tests;$MAKECMD clean_bin;time $MAKECMD;date) \
> $dir/gmake.log 2>&1
exit
else
echo "Unknwn option " $1
fi
fi
fi
fi
#
