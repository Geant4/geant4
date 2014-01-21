#!/bin/sh 

export DISPLAY_TMP=$DISPLAY
id=89
Njobs=0

# choose a free display number

while test $Njobs -ne 1; do
  (( id++ ))

  if test $id -gt 99 ; then
     echo "Cannot find a free display number, all /tmp/test202-[90-99].lock in use"
     exit 1
  fi

  /usr/bin/lockfile -l 3600 -r1 /tmp/test202-${id}.lock || continue

# check if there is a lock from Xserver for this display, if so try another
  [ -r /tmp/.X${id}-lock ] && continue
  pid=`cat /tmp/.X${id}-lock`
  ps -p $pid && continue

  XV_CMD="Xvfb :${id} -screen 0 1024x768x24 -nolisten tcp"
  $XV_CMD 2>&1 &
  XVFB_PID=$!
  sleep 1
  Njobs=`/bin/ps -f -p ${XVFB_PID} | /bin/grep -e " :${id} " | /usr/bin/wc -l `
  if  ( [ $Njobs -eq 1 ] ) ; then
# success quit search loop
     break
	  
  else  
     echo "Error: No Xserver for display :${id} started!"
     echo " is there a problem with lockfile/pipe still existing?"
	  echo "  Check for /tmp/.X${id}-lock and /tmp/.X11-unix/X${id}"
     rm -f /tmp/test202-${id}.lock
  fi	

done


export DISPLAY=:${id}
echo "launching $1"
$1 $2
rc=$?

kill  $XVFB_PID
rm -f /tmp/test202-${id}.lock

echo "end of test202 executable: RC=$rc"
exit $rc
