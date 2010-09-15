#!/bin/bash

# Program for running the same test case using both FORTRAN and C++
# versions. The script parses runs.conf file. The entries in the file
# are in format:

# name type A Z projectile energy events
# name = the name of the test (we use systematic runN format)
# type = full (INCL/ABLA), cascade (INCL only)
# A = mass number of the target 
# Z = charge number of the target 
# projectile = 1 (proton), 2 (neutron), 3 (pion+), 4 (pion-),
#              5(pion0), 6 (H2), 7 (H3), 8 (He3), 9 (He4)
# energy = projectile energy in MeV (100 - 3000 MeV)
# events = the number of events to be processed (includes transparent
#          events) 

# The output file names are constructed by this script and thus always
# use the same format. The format for the output files is:
# runN.root (C++ version)
# runN-ref.hbook (FORTRAN version)

# Environment variable $RUNFILE contains the name of the file where
# runs are defined. If this variable is NOT defined, the default run
# file is runs.conf.

if [ -n "$RUNFILE" ]; then
    echo "Using runfile: $RUNFILE"
    echo "RunID: $1"
else
    RUNFILE=runs.conf
fi

runID=$1

# Check that we have a valid runID
entries=$(cat $RUNFILE |grep "$runID\ " |wc -l)
if test $entries != "1"; then
    if test $entries = "0"; then
	echo -n "run.sh: Fatal error. Run "
	echo -n $runID
	echo " not defined."
    else
	echo -n "run.sh: Fatal error. More than one definition for run "
	echo $runID
	echo "Found definitions:"
	cat runs.conf |grep $runID
    fi
    exit 1
fi

name=$(cat $RUNFILE |grep "$runID\ " |awk '{print $1}')
type=$(cat $RUNFILE |grep "$runID\ " |awk '{print $2}') 
A=$(cat $RUNFILE |grep "$runID\ " |awk '{print $3}')
Z=$(cat $RUNFILE |grep "$runID\ " |awk '{print $4}')
proj=$(cat $RUNFILE |grep "$runID\ " |awk '{print $5}')
energy=$(cat $RUNFILE |grep "$runID\ " |awk '{print $6}')
events=$(cat $RUNFILE |grep "$runID\ " |awk '{print $7}')

datadir="tmp"
output=$datadir"/"$name".root"
refoutput=$datadir"/"$name"ref"
logfile=$datadir"/"$name".log"

# Run both FORTRAN and C++ tests by default
runcpp="1"
runfortran="1"
if test $# -ge 2; then
    if test $2 = "cpp"; then
	runcpp=1
	runfortran=0
    fi
    if test $2 = "fortran"; then
	runcpp=0
	runfortran=1
    fi
fi
if test $# -ge 3; then
    if test $3 = "cascade"; then
	type="cascade"
    fi
    if test $3 = "full"; then
	type="full"
    fi
fi

# Run the C++ standalone version:
if test $runcpp -eq 1; then
    echo "C++ standalone run started:" > $logfile
    date >> $logfile
    ./incl $type $A $Z $proj $energy $events $output >> $logfile
    echo "C++ standalone run complete:" >> $logfile
    date >> $logfile
fi

# Run the FORTRAN standalone version:
if test $runfortran -eq 1; then
    echo "---------------------------------" >> $logfile
    echo "FORTRAN standalone run started:" >> $logfile
    date >> $logfile
    ./runref.sh $type $A $Z $proj $energy $events $refoutput >> $logfile
    echo "FORTRAN standalone run complete:" >> $logfile
    date >> $logfile
fi
echo "---------------------------------" >> $logfile
echo "Done." >> $logfile
exit 0
