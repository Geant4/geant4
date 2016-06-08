#! /bin/csh -f

# check for LHC++ variables

if ( ! $?LHCXXTOP ) set LHCXXTOP = /afs/cern.ch/sw/lhcxx
if ( ! $?PLATF )    set PLATF    = $OS
if ( ! $?OSPACE_DIR ) \
   set OSPACE_DIR = $LHCXXTOP/specific/$PLATF/ObjectSpace/2.1/ToolKit

# check the LD_LIBRARY_PATH

if( $?LD_LIBRARY_PATH ) then
  (echo $LD_LIBRARY_PATH | grep HepODBMS >& /dev/null) || \
    setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:$HEP_ODBMS_DIR/lib
  if ( $PLATF == "Linux" ) then
    (echo $LD_LIBRARY_PATH | grep ObjectSpace >& /dev/null) || \
      setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:$OSPACE_DIR/lib
  endif
  echo LD_LIBRARY_PATH = $LD_LIBRARY_PATH
  echo ""
endif

