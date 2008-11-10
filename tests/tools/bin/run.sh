#!/bin/sh -m
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

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
    $G4STTDIR/bin/run.sh test01
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test01 1

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test02 2
    $G4STTDIR/bin/run.sh test02
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test02 2

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test02.hadron 3
  $G4STTDIR/bin/run.sh test02.hadron
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test02.hadron 3

###
###   Loop...
###
###    $G4STTDIR/bin/run.sh test05

#  ${G4STTDIR}/bin/geant4-unix.pl --start-test test06 4
#    $G4STTDIR/bin/run.sh test06
#  ${G4STTDIR}/bin/geant4-unix.pl --end-test test06 4

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test07 5
  $G4STTDIR/bin/run.sh test07
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test07 5

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test09 6
    $G4STTDIR/bin/run.sh test09
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test09 6

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test10 7
    $G4STTDIR/bin/run.sh test10
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test10 7

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test11 8
  $G4STTDIR/bin/run.sh test11
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test11 8

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test12 9
    $G4STTDIR/bin/run.sh test12
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test12 9

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test13 10
    $G4STTDIR/bin/run.sh test13
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test13 10

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test14 11
    $G4STTDIR/bin/run.sh test14
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test14 11

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test15 12
    $G4STTDIR/bin/run.sh test15
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test15 12

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test16 13
    $G4STTDIR/bin/run.sh test16
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test16 13

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test17 14
    $G4STTDIR/bin/run.sh test17
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test17 14

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test18 15
  $G4STTDIR/bin/run.sh test18
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test18 15

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test20 16
    $G4STTDIR/bin/run.sh test20
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test20 16

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test21 17
  $G4STTDIR/bin/run.sh test21
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test21 17

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test22 18
    $G4STTDIR/bin/run.sh test22
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test22 18

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test23 19
    $G4STTDIR/bin/run.sh test23
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test23 19

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test32 20
    $G4STTDIR/bin/run.sh test32
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test32 20

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test101 21
    $G4STTDIR/bin/run.sh test101
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test101 21

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test102 22
    $G4STTDIR/bin/run.sh test102
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test102 22

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test103 23
  $G4STTDIR/bin/run.sh test103
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test103 23

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test104 24
    $G4STTDIR/bin/run.sh test104
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test104 24

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test104.EMtest 25
    $G4STTDIR/bin/run.sh test104.EMtest
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test104.EMtest 25

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test105 26
    $G4STTDIR/bin/run.sh test105
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test105 26

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test106 27
    $G4STTDIR/bin/run.sh test106
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test106 27

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test501 28
    $G4STTDIR/bin/run.sh test501
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test501 28

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test502 29
    $G4STTDIR/bin/run.sh test502
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test502 29

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test503 30
    $G4STTDIR/bin/run.sh test503
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test503 30

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test504 31
    $G4STTDIR/bin/run.sh test504
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test504 31

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test505 32
    $G4STTDIR/bin/run.sh test505
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test505 32

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test506 33
    $G4STTDIR/bin/run.sh test506
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test506 33

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test507 34
    $G4STTDIR/bin/run.sh test507
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test507 34  

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test508 35
    $G4STTDIR/bin/run.sh test508
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test508 35
  
  ${G4STTDIR}/bin/geant4-unix.pl --start-test test509 36
    $G4STTDIR/bin/run.sh test509
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test509 36

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test5010 37
    $G4STTDIR/bin/run.sh test5010
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test5010 37  
  
  ${G4STTDIR}/bin/geant4-unix.pl --start-test test5011 38
    $G4STTDIR/bin/run.sh test5011
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test5011 38

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test5012 39
    $G4STTDIR/bin/run.sh test5012
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test5012 39 

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test5013 40
    $G4STTDIR/bin/run.sh test5013
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test5013 40

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test5014 41
    $G4STTDIR/bin/run.sh test5014
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test5014 41  


#  if [ $G4USE_HEPODBMS ] ; then
#   $G4STTDIR/bin/run.sh test401
#   $G4STTDIR/bin/run.sh test402
#  fi

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test601 42
    $G4STTDIR/bin/run.sh test601
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test601 42

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test602 43
    $G4STTDIR/bin/run.sh test602
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test602 43

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test701 44
    $G4STTDIR/bin/run.sh test701
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test701 44

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test702 45
    $G4STTDIR/bin/run.sh test702
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test702 45

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test703 46
    $G4STTDIR/bin/run.sh test703
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test703 46

###

#  ${G4STTDIR}/bin/geant4-unix.pl --start-test test801 41
#    $G4STTDIR/bin/run.sh test801
#  ${G4STTDIR}/bin/geant4-unix.pl --end-test test801 41
#
#  ${G4STTDIR}/bin/geant4-unix.pl --start-test test802 42
#    $G4STTDIR/bin/run.sh test802
#  ${G4STTDIR}/bin/geant4-unix.pl --end-test test802 42
#
#  ${G4STTDIR}/bin/geant4-unix.pl --start-test test803 43
#    $G4STTDIR/bin/run.sh test803
#  ${G4STTDIR}/bin/geant4-unix.pl --end-test test803 43
#
#  ${G4STTDIR}/bin/geant4-unix.pl --start-test test804 44
#    $G4STTDIR/bin/run.sh test804
#  ${G4STTDIR}/bin/geant4-unix.pl --end-test test804 44
#
#  ${G4STTDIR}/bin/geant4-unix.pl --start-test test805 45
#    $G4STTDIR/bin/run.sh test805
#  ${G4STTDIR}/bin/geant4-unix.pl --end-test test805 45
#
#  ${G4STTDIR}/bin/geant4-unix.pl --start-test test806 46
#    $G4STTDIR/bin/run.sh test806
#  ${G4STTDIR}/bin/geant4-unix.pl --end-test test806 46
#
#  ${G4STTDIR}/bin/geant4-unix.pl --start-test test807 47
#    $G4STTDIR/bin/run.sh test807
#  ${G4STTDIR}/bin/geant4-unix.pl --end-test test807 47


#  ${G4STTDIR}/bin/geant4-unix.pl --start-test test1001 48
#    $G4STTDIR/bin/run.sh test1001
#  ${G4STTDIR}/bin/geant4-unix.pl --end-test test1001 48

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test33 47
    $G4STTDIR/bin/run.sh test33
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test33 47
  
  ${G4STTDIR}/bin/geant4-unix.pl --start-test test33_1 48
    $G4STTDIR/bin/run.sh test33_1
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test33_1 48

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test34 49
    $G4STTDIR/bin/run.sh test34
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test34 49

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test24 50
    $G4STTDIR/bin/run.sh test24
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test24 50

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test25 51
    $G4STTDIR/bin/run.sh test25
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test25 51

#  ${G4STTDIR}/bin/geant4-unix.pl --start-test test26 52
#    $G4STTDIR/bin/run.sh test26
#  ${G4STTDIR}/bin/geant4-unix.pl --end-test test26 52


### Advanced examples


  ${G4STTDIR}/bin/geant4-unix.pl --start-test test107 53
    $G4STTDIR/bin/run.sh test107
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test107 53


  ${G4STTDIR}/bin/geant4-unix.pl --start-test test27 54
    $G4STTDIR/bin/run.sh test27
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test27 54


#############
  ${G4STTDIR}/bin/geant4-unix.pl --start-test test19 55
    $G4STTDIR/bin/run.sh test19
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test19 55

  ${G4STTDIR}/bin/geant4-unix.pl --start-test test29 56
    $G4STTDIR/bin/run.sh test29
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test29 56

  
  ${G4STTDIR}/bin/geant4-unix.pl --start-test test40 57
    $G4STTDIR/bin/run.sh test40
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test40 57
    
  ${G4STTDIR}/bin/geant4-unix.pl --start-test test28 58
    $G4STTDIR/bin/run.sh test28
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test28 58


  ${G4STTDIR}/bin/geant4-unix.pl --start-test test801 59
    $G4STTDIR/bin/run.sh test801
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test801 59
  
  ${G4STTDIR}/bin/geant4-unix.pl --start-test test500 60
    $G4STTDIR/bin/run.sh test500
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test500 60
  
  ${G4STTDIR}/bin/geant4-unix.pl --start-test test5015 61
  $G4STTDIR/bin/run.sh test5015
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test5015 61
  
  ${G4STTDIR}/bin/geant4-unix.pl --start-test test5016 62
  $G4STTDIR/bin/run.sh test5016
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test5016 62
  
  ${G4STTDIR}/bin/geant4-unix.pl --start-test test5017 63
  $G4STTDIR/bin/run.sh test5017
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test5017 63
  
  ${G4STTDIR}/bin/geant4-unix.pl --start-test test1002 64
  $G4STTDIR/bin/run.sh test1002
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test1002 64
  
  ${G4STTDIR}/bin/geant4-unix.pl --start-test test61 65
  $G4STTDIR/bin/run.sh test61
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test61 65
  
  ${G4STTDIR}/bin/geant4-unix.pl --start-test test6001 66
  $G4STTDIR/bin/run.sh test6001
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test6001 66
  
  ${G4STTDIR}/bin/geant4-unix.pl --start-test test05 67
  $G4STTDIR/bin/run.sh test05
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test05 67
  
  ${G4STTDIR}/bin/geant4-unix.pl --start-test test7001 68
  $G4STTDIR/bin/run.sh test7001
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test7001 68
  
  ${G4STTDIR}/bin/geant4-unix.pl --start-test test39 69
    $G4STTDIR/bin/run.sh test39
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test39 69
  
  ${G4STTDIR}/bin/geant4-unix.pl --start-test test2001 70
  $G4STTDIR/bin/run.sh test2001
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test2001 70
  
  ${G4STTDIR}/bin/geant4-unix.pl --start-test test2002 71
  $G4STTDIR/bin/run.sh test2002
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test2002 71
  
    ${G4STTDIR}/bin/geant4-unix.pl --start-test test2003 72
  $G4STTDIR/bin/run.sh test2003
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test2003 72
  
  ${G4STTDIR}/bin/geant4-unix.pl --start-test test5018 73
  $G4STTDIR/bin/run.sh test5018
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test5018 73
  
  ${G4STTDIR}/bin/geant4-unix.pl --start-test test704 74
  $G4STTDIR/bin/run.sh test704
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test704 74
  
  ${G4STTDIR}/bin/geant4-unix.pl --start-test test1004 75
  $G4STTDIR/bin/run.sh test1004
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test1004 75
  
  ${G4STTDIR}/bin/geant4-unix.pl --start-test test1008 76
  $G4STTDIR/bin/run.sh test1008
  ${G4STTDIR}/bin/geant4-unix.pl --end-test test1008 76
  
#  ${G4STTDIR}/bin/geant4-unix.pl --start-test test803 61
#    $G4STTDIR/bin/run.sh test803
#  ${G4STTDIR}/bin/geant4-unix.pl --end-test test803 61

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
        rm -f $dir/$1.exerciser$dot_G4LARGE_N.in
        $G4WORKDIR/bin/$G4SYSTEM/$shortname.hadronic.exerciser $G4LARGE_N \
        > $dir/$1.exerciser$dot_G4LARGE_N.in
        rm -f $dir/$1$dot_G4LARGE_N.out
        rm -f $dir/$1$dot_G4LARGE_N.err
        /usr/bin/time $G4WORKDIR/bin/$G4SYSTEM/$shortname \
        $dir/$1.exerciser$dot_G4LARGE_N.in \
        > $dir/$1$dot_G4LARGE_N.out 2> $dir/$1$dot_G4LARGE_N.err

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
# real RUN
#  ${G4STTDIR}/bin/geant4-unix.pl --start-test $shortname 1
          /usr/bin/time $G4WORKDIR/bin/$G4SYSTEM/$shortname \
          $G4INSTALL/tests/$shortname/$1$dot_G4LARGE_N.in \
          > $dir/$1$dot_G4LARGE_N.out 2> $dir/$1$dot_G4LARGE_N.err
#  ${G4STTDIR}/bin/geant4-unix.pl --end-test $shortname 1
#
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
