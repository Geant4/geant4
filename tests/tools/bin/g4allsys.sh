#!/usr/local/bin/bash
######################################
# Submit remote build&run 
######################################

############################################################
# 'account' MUST have next in PATH:
# /afs/cern.ch/sw/geant4/stt/ref[+]/src/geant4/tests/tools/bin 
# This is no loner required - but leave because there are side effects
#   to solve.
############################################################
account=stesting
 
while read host devprod debopt tag action actarg1 actarg2 actarg3 nonincremental
do
    sym_devprod=`echo $devprod|cut -c 1`
    if [ X$sym_devprod = Xd ]
    then
	REF=ref+
    else
	REF=ref
    fi

    sym_debopt=`echo $debopt|cut -c 1`
    if [ X$sym_debopt = Xo ]
    then
	g4sbr=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4/tests/tools/bin/g4sbr.o.sh
    else
	g4sbr=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4/tests/tools/bin/g4sbr.d.sh
    fi

    command="rsh -l $account $host $g4sbr $devprod $debopt $tag $action $actarg1 $actarg2 $actarg3 $nonincremental"
    if [ X$host = X\# ]
    then
	echo "String is comment!"
    else
	echo $command
	$command&
     fi	
done
