#! /bin/sh
# $Id: updt.sh,v 1.2 1999-01-20 10:28:28 allison Exp $
# Usage:
#   csh: updt.sh [-n] < something.sdb >& something.update.log
#    sh: updt.sh [-n] < something.sdb > something.update.log 2>&1

# .sdb files have a special format - see, e.g., geant4/tests/stt-ref.sdb.

while read module tag comments
do
#
# We have strange bug in CVS and NOT use -P flag!
#
    command="cvs $1 update -d  -r $tag $module"
    if [ $module = \# ]
    then
	echo $command - ignored
    else
	echo $command
	$command
    fi
done
