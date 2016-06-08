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

# ----- Objectivity variables setup

if ( ! $?OBJY_VERS ) then
  setenv OBJY_VERS 5.1
endif

if ( -r /afs/cern.ch/rd45/objectivity/objyenv.csh ) then 
  echo "Setting up OBJY_VERS $OBJY_VERS ..."
  source /afs/cern.ch/rd45/objectivity/objyenv.csh 
endif 

# ----- HepODBMS variables setup

setenv HEP_ODBMS_DIR       /afs/cern.ch/sw/lhcxx/specific/@sys/HepODBMS/0.3.0.1
setenv HEP_ODBMS_INCLUDES  ${HEP_ODBMS_DIR}/include
