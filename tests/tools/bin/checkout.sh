#! /bin/bash
# $Id$
# Usage: Checkout the HEAD of geant4/tests first.  Then from there:
#   csh: geant4/tests/tools/bin/checkout.sh \
           <beta01-01.sdb >&beta01-01.checkout.log
#    sh: geant4/tests/tools/bin/checkout.sh \
           <beta01-01.sdb >beta01-01.checkout.log 2>&1

while read module tag
do
    command="cvs co -r $tag $module"
    if [ $module = \# ]
    then
	echo $command - ignored
    else
	echo $command
	$command
    fi
done
