#!/usr/local/bin/bash
######################################
# Submit remote build&run 
######################################

############################################################
# 'account' MUST have next in PATH:
# /afs/cern.ch/sw/geant4/stt/prod[dev[1/2]]/src/geant4/tests/tools/bin 
# This is no longer required - but leave because there are side effects
#   to solve.
############################################################
account=stesting
 
while read host devprod debopt tag action actarg1 actarg2 actarg3 nonincremental
do

#   this full pathname yields a command line truncated by ps (spoiling g4status.sh's analysis)

    g4sbr=/afs/cern.ch/sw/geant4/stt/$devprod/src/geant4/tests/tools/bin/g4sbr.sh

    command="rsh -l $account $host $g4sbr $devprod $debopt $tag $action $actarg1 $actarg2 $actarg3 $nonincremental"
    if [ X$host = X ]
    then
        exit # end of requests on command line (is there any bash shell documentation?)
    fi
    if [ X$host = X\# ]
    then
	echo "String is comment: $host $devprod $debopt $tag $action $actarg1 $actarg2 $actarg3 $nonincremental"
    else
        echo $command
	echo $PATH >>  $host.$devprod.$debopt.log 2>&1 
	echo "______________________" >>  $host.$devprod.$debopt.log 2>&1 
	echo $command >>  $host.$devprod.$debopt.log 2>&1 
	$command >> $host.$devprod.$debopt.log 2>&1 &
    fi	
done
