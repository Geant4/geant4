#!/bin/sh -f
#!/usr/local/bin/bash 
#
#  Usage :
#     UNIX> cd <g4sttdir>/bin
#     UNIX> chmod u+x run.sh
#     UNIX> ./run.sh all
#     UNIX> ./run.sh test01
#
#  When env variables are correctly setted, this script could be 
# executed from everywhere. (That a bit anbitious, anywhere is already
#                            quite risky).
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

if [ -z "${G4STTDIR}" ] ; then
  echo "G4STTDIR not set. Execute first setup file !"
  exit
fi

##if [ -z "${LD_LIBRARY_PATH}" ] ; then
#    LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM
#    export LD_LIBARY_PATH
#    echo $LD_LIBRARY_PATH
##fi


# Make $G4WORKDIR/stt directory :
dir=$G4WORKDIR/stt
if [ ! -d $dir ] ; then
 echo "$dir does not exist."
 echo "This is now an error.  Create a directory or a symbolic link"
 echo "  to a directory by hand or in a calling script."
fi

# Make $G4WORKDIR/stt/$G4SYSTEM directory :
dir=$G4WORKDIR/stt/$G4SYSTEM
if [ ! -d $dir ] ; then
  mkdir $dir
  echo "$dir created."
fi

# Check setup

# I was wondering where these apparently useless directories were coming from.
#-- 
#-- if [ ! -d $G4INSTALL/tests/tools/lib/$G4SYSTEM ]; then
#--   mkdir $G4INSTALL/tests/tools/lib/$G4SYSTEM
#-- fi
#-- if [ ! -f $G4INSTALL/tests/tools/lib/$G4SYSTEM/libG4TestTool.a ]; then
#--   cd $G4INSTALL/tests/tools/lib; gmake
#-- fi
#-- if [ ! -d $G4INSTALL/tests/tools/bin/$G4SYSTEM ]; then
#--   mkdir $G4INSTALL/tests/tools/bin/$G4SYSTEM
#-- fi
#-- ##if [ ! -f $G4INSTALL/tests/tools/bin/$G4SYSTEM/diffhisto ]; then
#-- ##  cd $G4INSTALL/tests/tools/bin; gmake
#-- ##fi


if [ $1 = "all" ] ; then 

#
# Re-write it by looping against list
#

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test01 1
    nice $G4STTDIR/bin/run.sh test01
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test01 1

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test02 2
    nice $G4STTDIR/bin/run.sh test02
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test02 2

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test02.hadron 3
  nice $G4STTDIR/bin/run.sh test02.hadron
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test02.hadron 3

###
###   Loop...
###
###    nice $G4STTDIR/bin/run.sh test05

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test06 4
    nice $G4STTDIR/bin/run.sh test06
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test06 4

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test07 5
  nice $G4STTDIR/bin/run.sh test07
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test07 5

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test09 6
    nice $G4STTDIR/bin/run.sh test09
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test09 6

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test10 7
    nice $G4STTDIR/bin/run.sh test10
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test10 7
#
# NDL # Following tag 'neu-V03-02-02'
# On Thu, 26 Jul 2001, Hans-Peter Wellisch wrote:
# Please use G4NDL3.4 for test11, and G4NDL0.2 for the other tests.
# set default to G4NDL0.2 in setup scripts

  NeutronHPCrossSections=/afs/cern.ch/sw/geant4/dev/data/G4NDL3.5;export NeutronHPCrossSections
  echo "STT:NeutronHPCrossSections G4NDL3.5";

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test11 8
  nice $G4STTDIR/bin/run.sh test11
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test11 8

#  NeutronHPCrossSections=/afs/cern.ch/sw/geant4/dev/data/G4NDL0.2;export NeutronHPCrossSections
#  echo "STT:NeutronHPCrossSections G4NDL0.2";

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test12 9
    nice $G4STTDIR/bin/run.sh test12
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test12 9

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test13 10
    nice $G4STTDIR/bin/run.sh test13
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test13 10

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test14 11
    nice $G4STTDIR/bin/run.sh test14
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test14 11

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test15 12
    nice $G4STTDIR/bin/run.sh test15
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test15 12

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test16 13
    nice $G4STTDIR/bin/run.sh test16
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test16 13

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test17 14
    nice $G4STTDIR/bin/run.sh test17
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test17 14

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test18 15
  nice $G4STTDIR/bin/run.sh test18
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test18 15

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test20 16
    nice $G4STTDIR/bin/run.sh test20
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test20 16

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test21 17
  nice $G4STTDIR/bin/run.sh test21
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test21 17

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test22 18
    nice $G4STTDIR/bin/run.sh test22
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test22 18

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test23 19
    nice $G4STTDIR/bin/run.sh test23
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test23 19

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test32 20
    nice $G4STTDIR/bin/run.sh test32
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test32 20

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test101 21
    nice $G4STTDIR/bin/run.sh test101
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test101 21

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test102 22
    nice $G4STTDIR/bin/run.sh test102
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test102 22

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test103 23
  nice $G4STTDIR/bin/run.sh test103
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test103 23

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test104 24
    nice $G4STTDIR/bin/run.sh test104
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test104 24

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test104.EMtest 25
    nice $G4STTDIR/bin/run.sh test104.EMtest
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test104.EMtest 25

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test105 26
    nice $G4STTDIR/bin/run.sh test105
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test105 26

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test106 27
    nice $G4STTDIR/bin/run.sh test106
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test106 27

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test501 28
    nice $G4STTDIR/bin/run.sh test501
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test501 28

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test502 29
    nice $G4STTDIR/bin/run.sh test502
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test502 29

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test503 30
    nice $G4STTDIR/bin/run.sh test503
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test503 30

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test504 31
    nice $G4STTDIR/bin/run.sh test504
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test504 31

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test505 32
    nice $G4STTDIR/bin/run.sh test505
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test505 32

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test506 33
    nice $G4STTDIR/bin/run.sh test506
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test506 33

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test508 34
    nice $G4STTDIR/bin/run.sh test508
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test508 34

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test5010 35
    nice $G4STTDIR/bin/run.sh test5010
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test5010 35

  if [ $G4USE_HEPODBMS ] ; then
   nice $G4STTDIR/bin/run.sh test401
   nice $G4STTDIR/bin/run.sh test402
  fi

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test601 36
    nice $G4STTDIR/bin/run.sh test601
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test601 36

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test602 37
    nice $G4STTDIR/bin/run.sh test602
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test602 37

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test701 38
    nice $G4STTDIR/bin/run.sh test701
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test701 38

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test702 39
    nice $G4STTDIR/bin/run.sh test702
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test702 39

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test703 40
    nice $G4STTDIR/bin/run.sh test703
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test703 40

###

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test801 41
    nice $G4STTDIR/bin/run.sh test801
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test801 41

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test802 42
    nice $G4STTDIR/bin/run.sh test802
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test802 42

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test803 43
    nice $G4STTDIR/bin/run.sh test803
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test803 43

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test804 44
    nice $G4STTDIR/bin/run.sh test804
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test804 44

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test805 45
    nice $G4STTDIR/bin/run.sh test805
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test805 45

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test806 46
    nice $G4STTDIR/bin/run.sh test806
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test806 46

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test807 47
    nice $G4STTDIR/bin/run.sh test807
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test807 47

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test1001 48
    nice $G4STTDIR/bin/run.sh test1001
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test1001 48

else

#    echo n=$n
#    ${G4STTDIR}/bin/geant4-unix.pl --start-test $1 $n

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
    shortname=`basename $shortname .EMtest`


    timestamp=`date +%Y%m%d-%H%M%S`
    g4option=`echo $G4WORKDIR | sed "s!^.*${G4SYSTEM}/!!"`
    g4identify=`echo $G4WORKDIR | sed "s!^.*geant4/stt!stt!"`
    disabled=${G4INSTALL}/tests/disable.${shortname}.${G4SYSTEM}.${g4option}
    if [ -f $disabled ]; then
      echo "Disabled $1 in $g4identify $timestamp  $G4LARGE_N "
      echo "Disabled $1 in $G4WORKDIR $timestamp  $G4LARGE_N " > $dir/$1$dot_G4LARGE_N.err
      rm -f $dir/$1$dot_G4LARGE_N.out # or diff will run below
    elif [ ! -f $G4WORKDIR/bin/$G4SYSTEM/$shortname ]; then
      echo "Missing $1 in $g4identify $timestamp  $G4LARGE_N "
    else

      if [ $G4LARGE_N ]; then
        dot_G4LARGE_N=.$G4LARGE_N
      fi

      cd $G4INSTALL/tests/$shortname

      echo "Starting $1 in $g4identify $timestamp  $G4LARGE_N "

      if [ $1 = test02.hadron -o $1 = test11 -o $1 = test12 -o $1 = test13 \
        -o $1 = test15 -o $1 = test16 -o $1 = test21 ]
      then
        if [ $1 = test11 ]; then
          NeutronHPCrossSections=/afs/cern.ch/sw/geant4/dev/data/G4NDL3.5
          export NeutronHPCrossSections
          echo "STT:hadrons! NeutronHPCrossSections G4NDL3.5";
        fi
        rm -f $dir/$1.exerciser$dot_G4LARGE_N.in
        $G4WORKDIR/bin/$G4SYSTEM/$shortname.hadronic.exerciser $G4LARGE_N \
        > $dir/$1.exerciser$dot_G4LARGE_N.in
        rm -f $dir/$1$dot_G4LARGE_N.out
        rm -f $dir/$1$dot_G4LARGE_N.err
        /usr/bin/time $G4WORKDIR/bin/$G4SYSTEM/$shortname \
        $dir/$1.exerciser$dot_G4LARGE_N.in \
        > $dir/$1$dot_G4LARGE_N.out 2> $dir/$1$dot_G4LARGE_N.err
#        if [ $1 = test11 ]; then
#          NeutronHPCrossSections=/afs/cern.ch/sw/geant4/dev/data/G4NDL0.2
#          export NeutronHPCrossSections
#          echo "STT:hadrons! NeutronHPCrossSections G4NDL0.2";
#        fi

      elif [ $1 = test601 -o $1 = test602 ] 
      then
        if [ -z "$G4LARGE_N" -o \
          \( -n "$G4LARGE_N" -a \
             -f $G4INSTALL/tests/$shortname/$1$dot_G4LARGE_N.in \) ]; then
          rm -f $dir/$1$dot_G4LARGE_N.out
          rm -f $dir/$1$dot_G4LARGE_N.err
          /usr/bin/time $G4WORKDIR/bin/$G4SYSTEM/$shortname \
            $G4INSTALL/examples/extended/g3tog4/data/testmodel.dat \
            $G4INSTALL/tests/$shortname/$1$dot_G4LARGE_N.in \
          > $dir/$1$dot_G4LARGE_N.out 2> $dir/$1$dot_G4LARGE_N.err
        else
          echo "tests/$shortname/$1$dot_G4LARGE_N.in does not exist."
        fi

      else

        if [ -z "$G4LARGE_N" -o \
          \( -n "$G4LARGE_N" -a \
             -f $G4INSTALL/tests/$shortname/$1$dot_G4LARGE_N.in \) ]; then
          rm -f $dir/$1$dot_G4LARGE_N.out
          rm -f $dir/$1$dot_G4LARGE_N.err
          /usr/bin/time $G4WORKDIR/bin/$G4SYSTEM/$shortname \
          $G4INSTALL/tests/$shortname/$1$dot_G4LARGE_N.in \
          > $dir/$1$dot_G4LARGE_N.out 2> $dir/$1$dot_G4LARGE_N.err
        else
          echo "tests/$shortname/$1$dot_G4LARGE_N.in does not exist."
        fi
      fi

      timestamp=`date +%Y%m%d-%H%M%S`
      echo "Finished $1 in $g4identify $timestamp"

      if [ -f $1$dot_G4LARGE_N.evalsh ]; then
        rm -f $dir/$1$dot_G4LARGE_N.eval
        $1$dot_G4LARGE_N.evalsh $dir/$1$dot_G4LARGE_N.out \
        $1$dot_G4LARGE_N.out > $dir/$1$dot_G4LARGE_N.eval 2>&1
        if [ $? != 0 ]; then
          if [ ! -f $dir/$1$dot_G4LARGE_N.badflag ]; then
            touch $dir/$1$dot_G4LARGE_N.badflag
          fi
          echo $1$dot_G4LARGE_N run failure \
          >> $dir/$1$dot_G4LARGE_N.badflag 2>&1
        fi
      else
        rm -f $dir/$1$dot_G4LARGE_N.diff
        rm -f $dir/$1$dot_G4LARGE_N.diff_err
        if [ -f $dir/$1$dot_G4LARGE_N.out ]; then
          diff -w $dir/$1$dot_G4LARGE_N.out $1$dot_G4LARGE_N.out \
          > $dir/$1$dot_G4LARGE_N.diff 2> $dir/$1$dot_G4LARGE_N.diff_err
          #cat $dir/$1$dot_G4LARGE_N.diff
        fi
      fi

    fi

  fi

#  ${G4STTDIR}/bin/geant4-unix.pl --end-test $1 $n
#
#n=$[$n+1]
#echo n=$n

fi
