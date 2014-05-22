#!/bin/sh

# Test file name :
testFile1="visTest-DrawVolume_0.eps"
testFile2="visTest-SceneAdd_0.eps"
testFile3="visTest-SceneSet_0.eps"
testFile4="visTest-BeamOn_0.eps"

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
  echo "Check if output file exists ... "
  if ! test -f $output_folder/$testFile1; then
    echo "ERROR: $output_folder/$testFile1 not found. No output produced!"
    exit 1
  fi

  # mail me results
  mail -s "G4Testing visTest-DrawVolume_0.eps" garnier@lal.in2p3.fr < $output_folder/$testFile1
  mail -s "G4Testing visTest-VisSet_0.eps" garnier@lal.in2p3.fr < $output_folder/$testFile2
  mail -s "G4Testing visTest-VisAdd_0.eps" garnier@lal.in2p3.fr < $output_folder/$testFile3
  mail -s "G4Testing visTest-BeamOn_0.eps" garnier@lal.in2p3.fr < $output_folder/$testFile4
  sleep 10
  resCheck=0

for testFile in $testFile1 $testFile2 $testFile3 $testFile4
do
  if [`diff $reference_folder/$testFile $output_folder/$testFile | wc -l` == 0 ]; then
    echo "Check $output_folder/$testFile ....OK"
  else
    echo "ERROR: Output and reference are different for $output_folder/$testFile !"
    echo "======================================================"
    echo `diff $reference_folder/$testFile $output_folder/$testFile`
    echo "======================================================"
    mail -s "G4Testing ERROR Output and reference are different for $output_folder/$testFile" garnier@lal.in2p3.fr < $output_folder/$testFile
    resCheck=1
  fi

done

else
  echo "No GLX extension on Xvfb server"
  XVFB_PID="`pgrep -f "$XV_CMD"`"
  if [ XVFB_PID ]; then
    kill -9 $XVFB_PID
  fi
  exit 0
fi

exit $resCheck
