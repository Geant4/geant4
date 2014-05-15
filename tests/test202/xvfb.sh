#!/bin/sh

# Test file name :
testFile="visSceneTest_0.eps"

# reference folder :
reference_folder=$2-reference

# output folder : 
output_folder=$2-output

export DISPLAY_TMP=$DISPLAY
XV_CMD="Xvfb :99 -screen 0 1024x768x24 -nolisten tcp"

# is GLX extension present ?
rm -f XvfbExtensions
$XV_CMD 2>&1 &>XvfbExtensions &
sleep 1

# if GLX extension not present, return ok, but not run test
if grep "GLX" XvfbExtensions; then

  # OGLxx specific
  macro=$2.mac
  rm -rf $output_folder
  mkdir -p $output_folder

  # Launching
  echo "Launching application..."
  export DISPLAY=:99
  echo "launching $1"
  $1 $macro
  rc=$?

  export DISPLAY=$DISPLAY_TMP
  XVFB_PID="`pgrep -f "$XV_CMD"`"
  if [ $XVFB_PID ]; then
    kill -9 $XVFB_PID
  fi

  # check if all ok
  echo "Check if output and reference are similar... "
  if ! test -f $output_folder/$testFile; then
    echo "ERROR: $output_folder/$testFile not found. No output produced!"
    exit 1
  fi

  if [`diff $reference_folder/$testFile $output_folder/$testFile | wc -l` == 0 ]; then
    echo "All OK"
    exit 0
  else
    echo "ERROR: Output and reference are different!"
    echo "======================================================"
    echo `diff $reference_folder/$testFile $output_folder/$testFile`
    echo "======================================================"
    mail -s "ERROR: Output and reference are different for $2" garnier@lal.in2p3.fr < $output_folder/$testFile
    exit 1
  fi

else
  echo "No GLX extension on Xvfb server"
  XVFB_PID="`pgrep -f "$XV_CMD"`"
  if [ XVFB_PID ]; then
    kill -9 $XVFB_PID
  fi
  exit 0
fi
