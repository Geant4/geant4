#!/bin/sh -f
#
#
#
#  Usage :
#     UNIX> cd <g4install>/tests/tools/bin
#     UNIX> chmod u+x run.sh
#     UNIX> ./run.sh all
#     UNIX> ./run.sh test01
#
#  When env variables are correctly setted, this script could be 
# executed from everywhere.
#

# Some checks :
# set -x
if [ -z "${G4INSTALL}" ] ; then
  echo "G4INSTALL not set. Execute first setup file !"
  exit
fi

if [ -z "${G4WORKDIR}" ] ; then
  echo "G4WORKDIR not set. Execute first setup file !"
  exit
fi

if [ -z "${G4SYSTEM}" ] ; then
  echo "G4SYSTEM not set. Execute first setup file !"
  exit
fi

# Make $G4WORKDIR/stt directory :
dir=$G4WORKDIR/stt
if [ -d $dir ] ; then
  echo "" > /dev/null
else 
  mkdir $dir
  echo "$dir created."
fi

# Make $G4WORKDIR/stt/$G4SYSTEM directory :
dir=$G4WORKDIR/stt/$G4SYSTEM
if [ -d $dir ] ; then
  echo "" > /dev/null
else 
  mkdir $dir
  echo "$dir created."
fi

# Check setup
if [ ! -d $G4INSTALL/tests/tools/lib/$G4SYSTEM ]; then
  mkdir $G4INSTALL/tests/tools/lib/$G4SYSTEM
fi
if [ ! -f $G4INSTALL/tests/tools/lib/$G4SYSTEM/libG4TestTool.a ]; then
  cd $G4INSTALL/tests/tools/lib; gmake
fi
if [ ! -d $G4INSTALL/tests/tools/bin/$G4SYSTEM ]; then
  mkdir $G4INSTALL/tests/tools/bin/$G4SYSTEM
fi
##if [ ! -f $G4INSTALL/tests/tools/bin/$G4SYSTEM/diffhisto ]; then
##  cd $G4INSTALL/tests/tools/bin; gmake
##fi

if [ $1 = "all" ] ; then 

  nice $G4INSTALL/tests/tools/bin/run.sh test01
  nice $G4INSTALL/tests/tools/bin/run.sh test02
  nice $G4INSTALL/tests/tools/bin/run.sh test02.hadron
  nice $G4INSTALL/tests/tools/bin/run.sh test03
  nice $G4INSTALL/tests/tools/bin/run.sh test05
  nice $G4INSTALL/tests/tools/bin/run.sh test06
  nice $G4INSTALL/tests/tools/bin/run.sh test07
  nice $G4INSTALL/tests/tools/bin/run.sh test08
  nice $G4INSTALL/tests/tools/bin/run.sh test09
  nice $G4INSTALL/tests/tools/bin/run.sh test10
  nice $G4INSTALL/tests/tools/bin/run.sh test101
  nice $G4INSTALL/tests/tools/bin/run.sh test102
  nice $G4INSTALL/tests/tools/bin/run.sh test103
  nice $G4INSTALL/tests/tools/bin/run.sh test104
  nice $G4INSTALL/tests/tools/bin/run.sh test105
  nice $G4INSTALL/tests/tools/bin/run.sh test106
  nice $G4INSTALL/tests/tools/bin/run.sh test12
#  nice $G4INSTALL/tests/tools/bin/run.sh test13
# test13 tests kinetic model which is not yet available.
  nice $G4INSTALL/tests/tools/bin/run.sh test11
# test11 at end while it crashes on SUN in opt mode.

else

  if [ $1 = "test201" ] ; then 
  XENVIRONMENT=$G4INSTALL/tests/$1/$1.xrm;export XENVIRONMENT
    if [ $2 = "Xm" ] ; then 
      cd $G4INSTALL/tests/$1/basic
    else
      cd $G4INSTALL/tests/$1
    fi
    $G4WORKDIR/bin/$G4SYSTEM/$1 $2
    exit

  else

    # Other tests :
    shortname=`basename $1 .hadron`
    cd $G4INSTALL/tests/$shortname
    /bin/rm -f $dir/$1.out
    /bin/rm -f $dir/$1.diff
    if [ $1 = test02.hadron -o $1 = test11 -o $1 = test12 -o $1 = test13 ]
    then
      $G4WORKDIR/bin/$G4SYSTEM/$shortname.hadronic.exerciser \
      > $dir/$1.exerciser.in; \
      $G4WORKDIR/bin/$G4SYSTEM/$shortname \
      < $dir/$1.exerciser.in \
      > $dir/$1.out 2> $dir/$1.err
    else
      $G4WORKDIR/bin/$G4SYSTEM/$shortname \
      < $G4INSTALL/tests/$shortname/$1.in \
      > $dir/$1.out 2> $dir/$1.err
    fi
    diff $1.out $dir/$1.out > $dir/$1.diff 2> $dir/$1.diff_err
    #cat $dir/$1.diff

  fi

fi
