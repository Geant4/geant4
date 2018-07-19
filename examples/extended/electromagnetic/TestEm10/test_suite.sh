#!/bin/sh 

# run all TestEm10 macros and save their outputs in  log directory

CURDIR=`pwd`
OUTDIR=$CURDIR/log
PASSED="0"
FAILED="0"

if [ ! -d $OUTDIR ]; then
  mkdir -p $OUTDIR
fi

for MACRO in TestEm10.in alice06.mac bari05.mac  barr90.mac  harris73.mac salice.mac watase86.mac
do
  OUT=$OUTDIR/$MACRO.out

  echo  "... processing macro $MACRO"
  ./TestEm10 $MACRO > $OUT

  TMP_FAILED="0"
  if [ "$?" -ne "0" ]; then TMP_FAILED="1" ; fi
  if [ "$TMP_FAILED" -ne "0" ]; then FAILED=`expr $FAILED + 1`; else PASSED=`expr $PASSED + 1`; fi
done

# Print summary message
if [ "$FAILED" -eq "0" -a  "$PASSED" -ne "0" ]; then
  echo "... All ($PASSED) tests passed successfully."
elif [ "$FAILED" -ne "0" -a  "$PASSED" -eq "0" ]; then
  echo "... All ($FAILED) tests failed."
else
  echo "... $PASSED tests passed successfully."
  echo "... $FAILED tests failed."
fi
echo " "

cd $CURDIR

exit $FAILED
