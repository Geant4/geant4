#! /bin/sh

# check for LHC++ variables

if [ -z "$LHCXXTOP" ] ; then LHCXXTOP=/afs/cern.ch/sw/lhcxx ; fi
if [ -z "$PLATF" ]    ; then PLATF=$OS ; fi
if [ -z "$OSPACE_DIR" ] ; then
  OSPACE_DIR=$LHCXXTOP/specific/$PLATF/ObjectSpace/2.1/ToolKit
fi

# check the LD_LIBRARY_PATH

if [ -n "$LD_LIBRARY_PATH" ] ; then
  (echo $LD_LIBRARY_PATH | grep HepODBMS >& /dev/null) || \
   LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$HEP_ODBMS_DIR/lib; export LD_LIBRARY_PATH
  if [ $PLATF = "Linux" ] ; then
    (echo $LD_LIBRARY_PATH | grep ObjectSpace >& /dev/null) || \
     LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$OSPACE_DIR/lib; export LD_LIBRARY_PATH
  fi
  echo LD_LIBRARY_PATH = $LD_LIBRARY_PATH
  echo ""
fi

