#!/usr/bin/env bash
#
#############################################################################################################################
#####               This is an example of slurm file to  submit a job to execute dsbandrepair on cluster                #####
##### In this example, each node on cluster has memory of 31 Gb. And there are 16 cpus per node. In phycal stage, the   #####
##### parallel proceeses running on each node (buy setting num_ranks_pernode) is sett based the memory required to hold #####
##### geometry. For instance, with geometry set provided along this dsbandrepair, fbroblast and endothelium needs~4.5Gb,#####
#####  while yeast need ~4.5Gb and ~1.2GB respectively. Therefore, with a node of 31Gb memory, num_ranks_pernodeP       #####
##### can be set to 6 for fiborblast and endothelium, and 16 for yeast. For chemical stage, it is better to set the     #####    
##### number of chemical proceses on each node equal to the number of cpus. User should edit this file according to     #####
##### their need.                                                                                                       #####
#############################################################################################################################
##--------------------------------------------------------------------------------------------------------------------------##
#SBATCH --job-name="dsbandrepair"
#SBATCH --partition=std
#SBATCH --exclusive
#SBATCH --nodes=5  
##SBATCH --mem=190Gb 
##SBATCH --nodelist=node118
num_ranks_pernodeP=6 ## Number of Physical processes on each node. 
num_ranks_pernodeC=16  ## Number of chemical processes on all nodes. should be 1 cpu for 1 process
## You can change above setting for your need.
##--------------------------------------------------------------------------------------------------------------------------##
## for physical stage:
totalnumRankP=$(( $num_ranks_pernodeP*$SLURM_NNODES ))
## for chemicall stage:
totalnumRankC=$(( $num_ranks_pernodeC*$SLURM_NNODES ))
##--------------------------------------------------------------------------------------------------------------------------##
#" Check some requried files && folders
physmacfile="dsbandrepair.in"  #change it if you use other files for default
chemmacfile="chem.in" #change it if you use other files for default
flag="all"
##Read input arguments
for i in "$@"
do
    if [ $i = "-f" ] ;then shift;unset flag;flag=$1;shift;fi
    if [ $i = "-mP" ] ; then shift;unset physmacfile; physmacfile=$1;shift;fi
    if [ $i = "-mC" ] ; then shift;unset chemmacfile; chemmacfile=$1;shift;fi
done
##--------------------------------------------------------------------------------------------------------------------------##


logfolder="logs"
inputfolder="chem_input"
outputfolder="chem_output"
if [ ! -d $logfolder ]; then 
# folder to contain logfiles
    mkdir "$logfolder"
    mkdir "$logfolder/phys"
    mkdir "$logfolder/chem"
fi  
##--------------------------------------------------------------------------------------------------------------------------##
#START_TIME=$SECONDS
##--------------------------------------------------------------------------------------------------------------------------##
##If $flag = "phys", then run the physStage part
if [ $flag = "phys" ] || [ $flag = "all" ]; then          
    echo "Start running physical stage................."
    echo "This job will run with: "
    echo "=====> Number of nodes: $SLURM_NNODES"
    echo "=====> Number of Ranks: $numRanks"   
    mpiexec -np $totalnumRankP -npernode $num_ranks_pernodeP ./dsbandrepair $physmacfile 
    wait
    echo "End running physical stage................."
fi
##--------------------------------------------------------------------------------------------------------------------------##
##If $flag = "chem", then run the chemStage part
wait # make sure all above processes finish before chemStage starts
#sleep 1s
if [ $flag = "chem" ] || [ $flag = "all" ]; then
    if [ -d $logfolder/chem ]; then find $logfolder/chem/ -type f -delete;fi
    if [ -d $outputfolder ]; then
    find $outputfolder/ -type f -delete
    else
    mkdir $outputFolder
    fi
    echo "Start running chemical stage................."
    echo "with number of $totalnumRankC tasks.!!"
    echo "See $logfolder/* for running details"
    mpiexec -np $totalnumRankC ./dsbandrepair $chemmacfile chem $inputFolder > $logfolder/chem/log.dat
    wait
    echo "End running chemical stage on all tasks  ................."
fi
##--------------------------------------------------------------------------------------------------------------------------##
#echo "Elasped timed for $flag stage: $(($SECONDS - $START_TIME)) sec!!!"
##--------------------------------------------------------------------------------------------------------------------------##
