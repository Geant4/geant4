#!/bin/sh
# $Id: updt.sh,v 1.7 2002-07-21 19:16:07 stesting Exp $
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
	echo $module $tag $comments  - ignored
    else
	echo $command
	$command
    fi
done
