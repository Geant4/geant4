#! /bin/bash
# $Id: checkout.sh,v 1.1 1999-01-08 16:36:05 gunter Exp $
# Usage: Checkout the HEAD of geant4beta/tests first.  Then from there:
#   csh: geant4beta/tests/tools/bin/checkout.sh \
           <beta01-01.sdb >&beta01-01.checkout.log
#    sh: geant4beta/tests/tools/bin/checkout.sh \
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
