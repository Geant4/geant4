#
# G4ODBMS user example setup for CERN AFS
#
#  Put a line in your .cshrc or .login to source this file.
#
#  eg.  source ...<path of this file>.../g4odbms_setup.csh
#
#                                                   Youhei Morita
#

# ----- Objectivity variables setup

setenv OBJY_VERS 4.0.2

if ( -r /afs/cern.ch/rd45/objectivity/objyenv.csh ) then 
  source /afs/cern.ch/rd45/objectivity/objyenv.csh 
endif 

# ----- HepODBMS variables setup

setenv HEP_ODBMS_DIR        /afs/cern.ch/sw/lhcxx/specific/ibm/HepODBMS/0.0
setenv HEP_ODBMS_INCLUDES   ${HEP_ODBMS_DIR}/include
