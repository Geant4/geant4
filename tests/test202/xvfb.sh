#!/bin/sh
export DISPLAY_TMP=$DISPLAY
XV_CMD="Xvfb :99 -screen 0 1024x768x24 -nolisten tcp"
$XV_CMD 2>&1 &
export DISPLAY=:99
echo "launching $1"
$1 $2
rc=$?

export DISPLAY=$DISPLAY_TMP
echo "fin: RC=$rc"
XVFB_PID="`pgrep -f "$XV_CMD"`"
kill -9 $XVFB_PID
exit $rc
