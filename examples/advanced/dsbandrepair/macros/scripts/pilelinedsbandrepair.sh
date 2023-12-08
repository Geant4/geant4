#!/usr/bin/env bash
numRankP=1 # Number of ranks(cpu/cores) for physical stage ( 1 is default value)
numRankC=4 # Number of ranks(cores/cpus) for chemical stage( 4 is default value)
flag="all"; #default, 
inputFolder="chem_input"
outputFolder="chem_output"
physmacfile="dsbandrepair.in"  #change it if you use other files
chemmacfile="chem.in" #change it if you use other files
##--------------------------------------------------------------------------------------------------------------------------##
logfolder="logs"
if [ ! -d $logfolder ]; then 
# folder to contain logfiles
    mkdir "$logfolder"
    mkdir "$logfolder/phys"
    mkdir "$logfolder/chem"
fi  
##--------------------------------------------------------------------------------------------------------------------------##
##Read input arguments
for i in "$@"
do
    if [ $i = "-f" ] ;then shift;unset flag;flag=$1;shift;fi
    if [ $i = "-nRP" ] ;then shift;unset numRankP; numRankP=$1;shift;fi
    if [ $i = "-nRC" ] ; then shift;unset numRankC; numRankC=$1;shift;fi
    if [ $i = "-mP" ] ; then shift;unset physmacfile; physmacfile=$1;shift;fi
    if [ $i = "-mC" ] ; then shift;unset chemmacfile; chemmacfile=$1;shift;fi
done
##--------------------------------------------------------------------------------------------------------------------------##
echo "See $logfolder/* for running details"
#START_TIME=$SECONDS
##--------------------------------------------------------------------------------------------------------------------------##
#PhysStage:
if [ $flag = "all" ] || [ $flag = "phys" ]; then 
    if [ -d $logfolder/phys ]; then find $logfolder/phys/ -type f -delete;fi
    echo "Start running physical stage................."
    mpiexec -np $numRankP --bind-to none ./dsbandrepair $physmacfile > $logfolder/phys/log.dat
    wait
    echo "End running physical stage................."
fi
wait # make sure all above processes finish before chemStage starts

##--------------------------------------------------------------------------------------------------------------------------##

#ChemStage:
if [ $flag = "all" ] || [ $flag = "chem" ]; then
    if [ -d $logfolder/chem ]; then find $logfolder/chem/ -type f -delete;fi
    if [ -d $outputFolder ]; then
    find $outputFolder/ -type f -delete
    else
    mkdir "$outputFolder"
    fi
    echo "Start running chemical stage................."
    # Loop on each file of the $inputFolder folder
    mpiexec -np $numRankC --bind-to none ./dsbandrepair $chemmacfile chem $inputFolder > $logfolder/chem/log.dat
    wait
    echo "End running chemical stage................."
fi

##--------------------------------------------------------------------------------------------------------------------------##
#echo "Elasped timed for $flag stage: $(($SECONDS - $START_TIME)) sec!!!"
##--------------------------------------------------------------------------------------------------------------------------##