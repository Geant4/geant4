#! /bin/sh
# $Id: updt.sh,v 1.4 1999-07-28 10:37:23 stesting Exp $
# Usage:
#   csh: updt.sh < something.sdb >& something.update.log
#    sh: updt.sh < something.sdb > something.update.log 2>&1

# Uses environment variables NOTHING (-n) and DIRECTORIES (-d).

# .sdb files have a special format - see, e.g., geant4/tests/stt-prod.sdb.

while read module tag comments
do
#
# We have strange bug in CVS and NOT use -P flag!
#
    command="cvs $NOTHING update $DIRECTORIES -r $tag $module"
    if [ $module = \# ]
    then
	echo $command - ignored
    else
	echo $command
	$command
    fi
done
