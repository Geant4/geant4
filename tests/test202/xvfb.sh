#!/bin/sh
export DISPLAY_TMP=$DISPLAY
XV_CMD="Xvfb :99 -screen 0 1024x768x24 -nolisten tcp"
$XV_CMD&
export DISPLAY=:99
echo "launch"
$1 $2
export DISPLAY=$DISPLAY_TMP
echo "fin"
XVFB_PID="`pgrep -f "$XV_CMD"`"
kill -9 $XVFB_PID
exit 0
