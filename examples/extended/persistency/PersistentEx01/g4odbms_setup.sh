#
# G4ODBMS user example setup for CERN AFS
#
#  Put a line in your .zshrc to source this file.
#
#  eg.   . <path of this file>/g4odbms_setup.sh
#
#                                                   Youhei Morita
#

# ----- Use HepODBMS

G4USE_HEPODBMS=1; export G4USE_HEPODBMS

# ----- Objectivity variables setup

if [ -z "$OBJY_VERS" ] ; then
  OBJY_VERS=5.1
  export OBJY_VERS
fi

if [ -r /afs/cern.ch/rd45/objectivity/objyenv.sh ] ; then 
  echo "Setting up OBJY_VERS $OBJY_VERS ..."
  . /afs/cern.ch/rd45/objectivity/objyenv.sh 
fi 

# ----- HepODBMS variables setup

HEP_ODBMS_DIR=/afs/cern.ch/sw/lhcxx/specific/@sys/HepODBMS/0.3.0.1
export HEP_ODBMS_DIR
HEP_ODBMS_INCLUDES=${HEP_ODBMS_DIR}/include
export HEP_ODBMS_INCLUDES
