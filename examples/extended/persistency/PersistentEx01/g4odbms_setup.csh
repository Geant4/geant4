#
# G4ODBMS user example setup for CERN AFS
#
#  Put a line in your .cshrc or .login to source this file.
#
#  eg.  source ...<path of this file>.../g4odbms_setup.csh
#
#                                                   Youhei Morita
#

# ----- Use HepODBMS

setenv G4USE_HEPODBMS 1

# ----- LHCXX related variables

if( ! $?LHCXXTOP ) then
  setenv LHCXXTOP /afs/cern.ch/sw/lhcxx
endif

switch ($G4SYSTEM)
  case SUN-CC:
    setenv PLATF         sol7
    setenv CC_COMP       CC-5.2
    setenv OBJY_VERS     6.1.3
    setenv HEP_ODBMS_VER 0.3.3.1
    breaksw
  case Linux-g++:
    setenv PLATF         @sys
    setenv CC_COMP       gcc-2.95.2
    setenv OBJY_VERS     6.1.3
    setenv HEP_ODBMS_VER 0.3.3.1
    breaksw
  default:
    echo "Objectivity and HepODBMS versions not specified for G4SYSTEM = $G4SYSTEM."
    breaksw
endsw

if(! $?CC_COMP )then
  setenv CC_COMP ""
else
  if( "x${CC_COMP}" != "x" ) then
    setenv CC_COMP ${CC_COMP}/
  endif
endif

# ----- HepODBMS variables setup

setenv HEP_ODBMS_DIR $LHCXXTOP/specific/$PLATF/${CC_COMP}HepODBMS/${HEP_ODBMS_VER}
setenv HEP_ODBMS_INCLUDES  ${HEP_ODBMS_DIR}/include
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:$HEP_ODBMS_DIR/lib

# ----- Objectivity variables setup

if ( -r /afs/cern.ch/rd45/objectivity/objyenv.csh ) then 
  setenv OBJY_DIR $LHCXXTOP/specific/$PLATF/${CC_COMP}Objectivity/${OBJY_VERS}
  source /afs/cern.ch/rd45/objectivity/objyenv.csh 
else
  echo "/afs/cern.ch/rd45/objectivity/objyenv.csh not found!"
endif 

