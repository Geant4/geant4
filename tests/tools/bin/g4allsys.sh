#!/usr/local/bin/bash
######################################
# Submit remote build&run 
######################################

############################################################
# 'account' MUST have next in PATH:
# /afs/cern.ch/rd44/stt/ref+/src/geant4/tests/tools/bin 
############################################################
account=stesting
 
while read host flavour tag action actarg1 actarg2 actarg3 nonincremental
do
    sym_flavour=`echo $flavour|cut -c 1`
    if [ X$sym_flavour = Xo ]
    then
	g4sbr=g4sbr.o.sh
    else
	g4sbr=g4sbr.d.sh
    fi
    command="rsh -l $account $host $g4sbr $sym_flavour $tag $action $actarg1 $actarg2 $actarg3 $nonincremental"
    if [ X$host = X\# ]
    then
	echo "String is comment!"
    else
	echo $command
	$command&
     fi	
done
