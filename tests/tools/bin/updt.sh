#! /bin/sh
# $Id: updt.sh,v 1.6 2000-01-17 09:46:45 stesting Exp $
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
# We are now at CVS 1.10.3.  -P re-instated!
#
    command="cvs $NOTHING update -d -P -r $tag $module"
    if [ $module = \# ]
    then
	echo $command - ignored
    else
	echo $command
	$command
    fi
done
