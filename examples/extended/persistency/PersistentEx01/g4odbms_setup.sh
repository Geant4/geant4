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

export G4USE_HEPODBMS=1

# ----- LHCXX related variables

if [ -z "$LHCXXTOP" ] ; then
  export LHCXXTOP=/afs/cern.ch/sw/lhcxx
fi

export CC_COMP=""
case $G4SYSTEM in
SUN-CC5)
    export PLATF=sol7
    export CC_COMP=CC-5.2/
    export OBJY_VERS=6.1.3
    export HEP_ODBMS_VER=0.3.2.10
    ;;
Linux-g++)
    export PLATF=@sys
    export CC_COMP=gcc-2.95.2/
    export OBJY_VERS=6.0
    export HEP_ODBMS_VER=0.3.2.3
    ;;
default)
    echo "Objectivity and HepODBMS versions not specified for G4SYSTEM = $G4SYSTEM."
    ;;
esac

# ----- HepODBMS variables setup

export HEP_ODBMS_DIR=$LHCXXTOP/specific/$PLATF/${CC_COMP}HepODBMS/${HEP_ODBMS_VER}
export HEP_ODBMS_INCLUDES=${HEP_ODBMS_DIR}/include
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$HEP_ODBMS_DIR/lib

# ----- Objectivity variables setup

if [ -r /afs/cern.ch/rd45/objectivity/objyenv.sh ] ; then 
  export OBJY_DIR=$LHCXXTOP/specific/$PLATF/${CC_COMP}Objectivity/${OBJY_VERS}
  . /afs/cern.ch/rd45/objectivity/objyenv.sh
fi

