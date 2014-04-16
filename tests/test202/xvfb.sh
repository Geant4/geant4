#!/bin/sh
export DISPLAY_TMP=$DISPLAY
XV_CMD="Xvfb :99 -screen 0 1024x768x24 -nolisten tcp"
# is GLX extension present ?
rm -f XvfbExtensions
$XV_CMD 2>&1 &>XvfbExtensions &
sleep 1
# if GLX extension not present, return ok, but not run test
if grep "GLX" XvfbExtensions; then
  export DISPLAY=:99
  echo "launching $1"
  $1 $2
  rc=$?

  export DISPLAY=$DISPLAY_TMP
  XVFB_PID="`pgrep -f "$XV_CMD"`"
  if [ $XVFB_PID ]; then
    kill -9 $XVFB_PID
  fi
  exit $rc
else
  echo "No GLX extension on Xvfb server"
  XVFB_PID="`pgrep -f "$XV_CMD"`"
  if [ XVFB_PID ]; then
    kill -9 $XVFB_PID
  fi
  exit 0
fi
