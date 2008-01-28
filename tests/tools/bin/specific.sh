#
# This file is sourced by setup.sh script.
# 
####################################################
####################################################
### Specific STT member machines configurations  ###
####################################################
####################################################
#
# Specfy CLHEP 1.5 as new, run persistance tests 401,402
# Use /tmp for G4TMP on dxplus01,02,03,04,
# Use /tmp for G4TMP on all dxplus, retain dev1/dev2/prod identity
#

REF=undefined
if [ `pwd | grep /stt/dev1/` ]; then
export  REF=dev1
fi
if [ `pwd | grep /stt/dev2/` ]; then
export  REF=dev2
fi
if [ `pwd | grep /stt/prod/` ]; then
export  REF=prod
fi

#export CLHEP_VERSION=pro
export CLHEP_VERSION=2.0.3.2

if [ $G4DEBUG ]; then
export DEBOPT=debug
else
export DEBOPT=optim
fi

# General G4 build flags :
export G4UI_BUILD_TERMINAL_SESSION=1
export G4UI_BUILD_GAG_SESSION=1
export G4VIS_BUILD_DAWN_DRIVER=1
export G4VIS_BUILD_DAWNFILE_DRIVER=1
export G4VIS_BUILD_VRML_DRIVER=1
export G4VIS_BUILD_VRMLFILE_DRIVER=1
export G4VIS_BUILD_RAYTRACER_DRIVER=1
export G4VIS_BUILD_RAYTRACERX_DRIVER=1
export G4VIS_BUILD_ASCIITREE_DRIVER=1

export G4LIB_BUILD_G3TOG4=1
export G4USE_G3TOG4=1

export G4LIB_BUILD_GDML=1
export G4LIB_USE_GDML=1


UNAMEN=`uname -n `
echo UNAMEN $UNAMEN

if [ `uname -n | grep sunasd1` ]; then
  export G4USE_OSPACE=1
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
  export G4SYSTEM=SUN-CC
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4STTDIR=/afs/cern.ch/sw/geant4/stt/$REF/testtools/geant4/tests/tools
  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4LIB=$G4WORKDIR/lib
  # G4 build flags :
  #######export G4UI_BUILD_XM_SESSION=1
  #######export G4VIS_BUILD_OPENGLXM_DRIVER=1
  #######export G4VIS_BUILD_OPENGLX_DRIVER=1
  #######export G4VIS_BUILD_OIX_DRIVER=1
fi

if [ `uname -n | grep refsol7` ]; then
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4STTDIR=/afs/cern.ch/sw/geant4/stt/$REF/testtools/geant4/tests/tools
  if [ $G4STTNONISO ]; then
    echo "refsol7 only supports the ISO compiler"
    export G4SYSTEM=SUN-CC
    export DEBOPT=${DEBOPT}_NONISO
    export G4USE_OSPACE=1
    export PATH=`echo $PATH | sed s/SUNWspro50/SUNWspro/`
    export CLHEP_BASE_DIR=/afs/cern.ch/sw/geant4/dev/CLHEP/SUN-CC/$CLHEP_VERSION
    # Persistency...
    if [ X$G4USE_HEPODBMS = Xaaa ]; then
       . $G4INSTALL/examples/extended/persistency/PersistentEx01/g4odbms_setup.sh
       export G4EXAMPLE_FDID=207
    fi
  else
# G4SYSTEM was changed!
    export G4SYSTEM=SUN-CC
    export DEBOPT=${DEBOPT}_53
    unset G4USE_OSPACE
    export PATH=/opt/forte_6.2/SUNWspro/bin/:${PATH}
#    export PATH=`echo $PATH | sed s/SUNWspro50/SUNWspro/`
    export LD_LIBRARY_PATH=/opt/forte_6.2/SUNWspro/lib/:${LD_LIBRARY_PATH}
    # No Persistency tonight ...
    unset G4USE_HEPODBMS
    export CLHEP_BASE_DIR=/afs/cern.ch/sw/geant4/dev/CLHEP/SUN-CC/$CLHEP_VERSION
  fi
  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4LIB=$G4WORKDIR/lib
  # G4 build flags :
  #######export G4UI_BUILD_XM_SESSION=1
  #######export G4VIS_BUILD_OPENGLXM_DRIVER=1
  export G4VIS_BUILD_OPENGLX_DRIVER=1
  export OGLHOME=/usr/local
  export OGLFLAGS="-I$OGLHOME/include"
  export OGLLIBS="-L$OGLHOME/lib -lMesaGLU -lMesaGL"
  #######export G4VIS_BUILD_OIX_DRIVER=1
fi

if [ `uname -n | grep refsol8` ]; then
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4STTDIR=/afs/cern.ch/sw/geant4/stt/$REF/testtools/geant4/tests/tools
    export G4SYSTEM=SUN-CC
    export DEBOPT=${DEBOPT}_54

    unset G4USE_OSPACE

    export PATH=/afs/cern.ch/project/sun/solaris/opt/SUNWspro8/bin:${PATH}
    export PATH=/usr/local/gcc-alt-3.2.3/bin:${PATH}
#    export PATH=`echo $PATH | sed s/SUNWspro50/SUNWspro/`
    export LD_LIBRARY_PATH=/afs/cern.ch/project/sun/solaris/opt/SUNWspro8/lib:${LD_LIBRARY_PATH}
    export LD_LIBRARY_PATH=/usr/local/gcc-alt-3.2.3/lib:${LD_LIBRARY_PATH}
    # No Persistency tonight ...
    unset G4USE_HEPODBMS
#    export CLHEP_BASE_DIR=/afs/cern.ch/sw/geant4/dev/CLHEP/pro/SUN-CC
#    export CLHEP_BASE_DIR=/afs/cern.ch/sw/geant4/dev/CLHEP/1.9.1.2/sunos58_CC54

  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4LIB=$G4WORKDIR/lib
# Take CLHEP with links to lcg area
  export CLHEP_BASE_DIR=$G4WORKDIR/clhep

  # G4 build flags :
  #######export G4UI_BUILD_XM_SESSION=1
  #######export G4VIS_BUILD_OPENGLXM_DRIVER=1
  export G4VIS_BUILD_OPENGLX_DRIVER=1
#  export G4VIS_USE_OPENGLX=1
  export OGLHOME=/usr/local
  export OGLFLAGS="-I$OGLHOME/include"
  export OGLLIBS="-L$OGLHOME/lib -lMesaGLU -lMesaGL"
  #######export G4VIS_BUILD_OIX_DRIVER=1
fi

if [ `uname -n | grep refsol9` ]; then
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4STTDIR=/afs/cern.ch/sw/geant4/stt/$REF/testtools/geant4/tests/tools
    export G4SYSTEM=SUN-CC
    export DEBOPT=${DEBOPT}_55

    unset G4USE_OSPACE

    export PATH=/afs/cern.ch/project/sun/solaris/opt/SUNWspro8/bin:${PATH}
    export PATH=/usr/local/gcc-alt-3.2.3/bin:${PATH}
#    export LD_LIBRARY_PATH=/afs/cern.ch/project/sun/solaris/opt/SUNWspro8/lib:${LD_LIBRARY_PATH}
#    export LD_LIBRARY_PATH=/usr/local/gcc-alt-3.2.3/lib:${LD_LIBRARY_PATH}
    # No Persistency tonight ...
    unset G4USE_HEPODBMS

  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4LIB=$G4WORKDIR/lib
# Take CLHEP with links to lcg area
  export CLHEP_BASE_DIR=$G4WORKDIR/clhep

  # G4 build flags :
  #######export G4UI_BUILD_XM_SESSION=1
  #######export G4VIS_BUILD_OPENGLXM_DRIVER=1
  export G4VIS_BUILD_OPENGLX_DRIVER=1
#  export G4VIS_USE_OPENGLX=1
  export OGLHOME=/usr/local
  export OGLFLAGS="-I$OGLHOME/include"
  export OGLLIBS="-L$OGLHOME/lib -lMesaGLU -lMesaGL"
  #######export G4VIS_BUILD_OIX_DRIVER=1
fi

if [ `uname -n | grep tersk08` ]; then
  export CVSROOT=:ext:stesting@Geant4.cvs.cern.ch:/cvs/Geant4
  export CVS_RSH=ssh
  export G4INSTALL=/afs/slac.stanford.edu/u/sf/stesting/stt/$REF/src/geant4
  export G4STTDIR=/afs/slac.stanford.edu/u/sf/stesting/stt/$REF/testtools/geant4/tests/tools
    export G4SYSTEM=SUN-CC
    export DEBOPT=${DEBOPT}_59_58

    unset G4USE_OSPACE

    # No Persistency tonight ...
    unset G4USE_HEPODBMS
#    export CLHEP_BASE_DIR=/afs/cern.ch/sw/geant4/dev/CLHEP/SUN-CC/$CLHEP_VERSION

  export G4WORKDIR=/afs/slac.stanford.edu/u/sf/stesting/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4LIB=$G4WORKDIR/lib
# Take CLHEP with links to lcg area
  export CLHEP_BASE_DIR=$G4WORKDIR/clhep
  
  # G4 build flags :
  #######export G4UI_BUILD_XM_SESSION=1
  #######export G4VIS_BUILD_OPENGLXM_DRIVER=1
  export G4VIS_BUILD_OPENGLX_DRIVER=1
#  export G4VIS_USE_OPENGLX=1
  export OGLHOME=/usr/local
  export OGLFLAGS="-I$OGLHOME/include"
  export OGLLIBS="-L$OGLHOME/lib -lMesaGLU -lMesaGL"
  #######export G4VIS_BUILD_OIX_DRIVER=1
  
  export G4LEDATA=$G4WORKDIR/data/G4EMLOW
  export G4NEUTRONHPDATA=$G4WORKDIR/data/G4NDL
  export NeutronHPCrossSections=$G4WORKDIR/data/G4NDL
  export G4LEVELGAMMADATA=$G4WORKDIR/data/PhotonEvaporation
  export G4ELASTICDATA=$G4WORKDIR/data/G4ELASTIC
  export G4RADIOACTIVEDATA=$G4WORKDIR/data/RadioactiveDecay
  
fi


if [ `uname -n | grep refsol7AAA` ]; then
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4STTDIR=/afs/cern.ch/sw/geant4/stt/$REF/testtools/geant4/tests/tools
  if [ $G4STTNONISO ]; then
    echo "refsol7 only supports the ISO compiler"
    export G4SYSTEM=SUN-CC
    export DEBOPT=${DEBOPT}_NONISO
    export G4USE_OSPACE=1
    export PATH=`echo $PATH | sed s/SUNWspro50/SUNWspro/`
    export CLHEP_BASE_DIR=/afs/cern.ch/sw/geant4/dev/CLHEP/SUN-CC/$CLHEP_VERSION
    # Persistency...
    if [ X$G4USE_HEPODBMS = Xaaa ]; then
       . $G4INSTALL/examples/extended/persistency/PersistentEx01/g4odbms_setup.sh
       export G4EXAMPLE_FDID=207
    fi
  else
# G4SYSTEM was changed!
    export G4SYSTEM=SUN-CC
    export DEBOPT=${DEBOPT}_ISO
    unset G4USE_OSPACE
    export PATH=`echo $PATH | sed s/SUNWspro50/SUNWspro/`
    # No Persistency tonight ...
    unset G4USE_HEPODBMS
    export CLHEP_BASE_DIR=/afs/cern.ch/sw/geant4/dev/CLHEP/SUN-CC/$CLHEP_VERSION
  fi
  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4LIB=$G4WORKDIR/lib
  # G4 build flags :
  #######export G4UI_BUILD_XM_SESSION=1
  #######export G4VIS_BUILD_OPENGLXM_DRIVER=1
  export G4VIS_BUILD_OPENGLX_DRIVER=1
  export OGLHOME=/usr/local
  export OGLFLAGS="-I$OGLHOME/include"
  export OGLLIBS="-L$OGLHOME/lib -lMesaGLU -lMesaGL"
  #######export G4VIS_BUILD_OIX_DRIVER=1
fi

# BLOCKed!
if [ `uname -n | grep sundev` ]; then
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4STTDIR=/afs/cern.ch/sw/geant4/stt/$REF/testtools/geant4/tests/tools
  if [ $G4STTNONISO ]; then
    echo "refsol7 only supports the ISO compiler"
    export G4SYSTEM=SUN-CC
    export DEBOPT=${DEBOPT}_NONISO
    export G4USE_OSPACE=1
    export PATH=`echo $PATH | sed s/SUNWspro50/SUNWspro/`
    export CLHEP_BASE_DIR=/afs/cern.ch/sw/geant4/dev/CLHEP/SUN-CC/$CLHEP_VERSION
    # Persistency...
    if [ X$G4USE_HEPODBMS = Xaaa ]; then
       . $G4INSTALL/examples/extended/persistency/PersistentEx01/g4odbms_setup.sh
       export G4EXAMPLE_FDID=207
    fi
  else
    export G4SYSTEM=SUN-CC
    export DEBOPT=${DEBOPT}_ISO
    unset G4USE_OSPACE
    export PATH=`echo $PATH | sed s/SUNWspro50/SUNWspro/`
    # No Persistency tonight ...
    unset G4USE_HEPODBMS
    export CLHEP_BASE_DIR=/afs/cern.ch/sw/geant4/dev/CLHEP/SUN-CC/$CLHEP_VERSION
  fi
  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4LIB=$G4WORKDIR/lib
  # G4 build flags :
  #######export G4UI_BUILD_XM_SESSION=1
  #######export G4VIS_BUILD_OPENGLXM_DRIVER=1
  export G4VIS_BUILD_OPENGLX_DRIVER=1
  export OGLHOME=/usr/local
  export OGLFLAGS="-I$OGLHOME/include"
  export OGLLIBS="-L$OGLHOME/lib -lMesaGLU -lMesaGL"
  #######export G4VIS_BUILD_OIX_DRIVER=1
fi



if [ `uname -n | grep suncmsb` ]; then
  export G4USE_OSPACE=1
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
  export G4SYSTEM=SUN-CC
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4STTDIR=/afs/cern.ch/sw/geant4/stt/$REF/testtools/geant4/tests/tools
  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4LIB=$G4WORKDIR/lib
#
  export CLHEP_BASE_DIR=/afs/cern.ch/user/s/stesting/work/clhep
  export CLHEP_LIB=CLHEP-CC
  export RWBASE=/afs/cern.ch/user/s/stesting/work/rogue
  # G4 build flags :
  #######export G4UI_BUILD_XM_SESSION=1
  #######export G4VIS_BUILD_OPENGLXM_DRIVER=1
  #######export G4VIS_BUILD_OPENGLX_DRIVER=1
  #######export G4VIS_BUILD_OIX_DRIVER=1
fi

if [ `uname -n | grep sungeant` ]; then
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4STTDIR=/afs/cern.ch/sw/geant4/stt/$REF/testtools/geant4/tests/tools
  if [ $G4STTNONISO ]; then
    export G4SYSTEM=SUN-CC4
    export DEBOPT=${DEBOPT}_NONISO
    export G4USE_OSPACE=1
    export PATH=`echo $PATH | sed s/SUNWspro50/SUNWspro/`
###TEMPORARY!
    unset G4USE_HEPODBMS
#    export CLHEP_BASE_DIR=/afs/cern.ch/sw/geant4/dev/CLHEP/SUN-CC4/$CLHEP_VERSION
    # Persistency...
    if [ X$G4USE_HEPODBMS = Xaaa ]; then
       . $G4INSTALL/examples/extended/persistency/PersistentEx01/g4odbms_setup.sh
       export G4EXAMPLE_FDID=207
    fi
  else
    export G4SYSTEM=SUN-CC
    export DEBOPT=${DEBOPT}_ISO
    unset G4USE_OSPACE
    . /afs/cern.ch/sw/geant4/dev/scripts/CC8.sh
#    export PATH=/afs/cern.ch/project/sun/solaris/opt/SUNWspro7/bin:${PATH}
#    export PATH=`echo $PATH | sed s/SUNWspro50/SUNWspro/`
#   export LD_LIBRARY_PATH=/afs/cern.ch/project/sun/solaris/opt/SUNWspro7/lib:${LD_LIBRARY_PATH}
#    export PATH=`echo $PATH | sed s/SUNWspro/SUNWspro50/`
#    export PATH=`echo $PATH | sed s/SUNWspro50/SUNWspro/`
    # No Persistency...
    unset G4USE_HEPODBMS
    # Persistency...
#   if [ X$G4USE_HEPODBMS = Xaaa ]; then
    if [ X$G4USE_HEPODBMS = Xaaa ]; then
       . $G4INSTALL/examples/extended/persistency/PersistentEx01/g4odbms_setup.sh
       export G4EXAMPLE_FDID=207
    fi
#    export CLHEP_BASE_DIR=/afs/cern.ch/sw/geant4/dev/CLHEP/SUN-CC/$CLHEP_VERSION
#    export CLHEP_BASE_DIR=/afs/cern.ch/sw/geant4/dev/CLHEP/1.9.1.2/SUN-CC
  fi
  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4LIB=$G4WORKDIR/lib
# Take CLHEP with links to lcg area
  export CLHEP_BASE_DIR=$G4WORKDIR/clhep
  
  export XERCESCROOT=/afs/cern.ch/sw/geant4/dev/XercesC/2.8.0/sunos58_CC54
  
  export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$XERCESCROOT/lib
  
  # G4 build flags :
  #######export G4UI_BUILD_XM_SESSION=1
  #######export G4VIS_BUILD_OPENGLXM_DRIVER=1
  export G4VIS_BUILD_OPENGLX_DRIVER=1
  export OGLHOME=/usr/local
  export OGLFLAGS="-I$OGLHOME/include"
  export OGLLIBS="-L$OGLHOME/lib -lMesaGLU -lMesaGL"
  #######export G4VIS_BUILD_OIX_DRIVER=1
  export G4VIS_BUILD_OPENGLX_DRIVER=1
#  export G4VIS_USE_OPENGLX=1

fi

# refsol8 -> sundev008
if [ `uname -n | grep sundev008` ]; then
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4STTDIR=/afs/cern.ch/sw/geant4/stt/$REF/testtools/geant4/tests/tools
  if [ $G4STTNONISO ]; then
    export G4SYSTEM=SUN-CC4
    export DEBOPT=${DEBOPT}_NONISO
    export G4USE_OSPACE=1
    export PATH=`echo $PATH | sed s/SUNWspro50/SUNWspro/`
###TEMPORARY!
    unset G4USE_HEPODBMS
    export CLHEP_BASE_DIR=/afs/cern.ch/sw/geant4/dev/CLHEP/SUN-CC4/$CLHEP_VERSION
    # Persistency...
    if [ X$G4USE_HEPODBMS = Xaaa ]; then
       . $G4INSTALL/examples/extended/persistency/PersistentEx01/g4odbms_setup.sh
       export G4EXAMPLE_FDID=207
    fi
  else
    export G4SYSTEM=SUN-CC
    export DEBOPT=${DEBOPT}_ISO
    unset G4USE_OSPACE
    export PATH=`echo $PATH | sed s/SUNWspro/SUNWspro50/`
    export PATH=`echo $PATH | sed s/SUNWspro50/SUNWspro/`
    # No Persistency...
    unset G4USE_HEPODBMS
    # Persistency...
#   if [ X$G4USE_HEPODBMS = Xaaa ]; then
    if [ X$G4USE_HEPODBMS = Xaaa ]; then
       . $G4INSTALL/examples/extended/persistency/PersistentEx01/g4odbms_setup.sh
       export G4EXAMPLE_FDID=207
    fi
#    export CLHEP_BASE_DIR=/afs/cern.ch/sw/geant4/dev/CLHEP/SUN-CC/$CLHEP_VERSION
  fi
  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4LIB=$G4WORKDIR/lib
# Take CLHEP with links to lcg area
  export CLHEP_BASE_DIR=$G4WORKDIR/clhep  
  # G4 build flags :
  #######export G4UI_BUILD_XM_SESSION=1
  #######export G4VIS_BUILD_OPENGLXM_DRIVER=1
  export G4VIS_BUILD_OPENGLX_DRIVER=1
  export OGLHOME=/usr/local
  export OGLFLAGS="-I$OGLHOME/include"
  export OGLLIBS="-L$OGLHOME/lib -lMesaGLU -lMesaGL"
  #######export G4VIS_BUILD_OIX_DRIVER=1
fi


if [ $UNAMEN = pcgeant ]; then
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
  export G4SYSTEM=Linux-egcs
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4STTDIR=/afs/cern.ch/sw/geant4/stt/$REF/testtools/geant4/tests/tools
  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4LIB=$G4WORKDIR/lib
  export CLHEP_BASE_DIR=/afs/cern.ch/sw/geant4/dev/CLHEP/Linux-egcs/$CLHEP_VERSION
#
# The G4SYSTEm (and directory) was changed!
#  export CLHEP_LIB=CLHEP-egcs

  export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${CLHEP_BASE_DIR}/lib:${LD_LIBRARY_PATH}
  
  # G4 build flags :
  ######export G4UI_BUILD_XM_SESSION=1
#  export G4VIS_BUILD_OPENGLXM_DRIVER=1
#  export G4VIS_BUILD_OPENGLX_DRIVER=1
#  export G4VIS_BUILD_OIX_DRIVER=1
fi

if [ $UNAMEN = pcitapi08 ]; then
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
  export G4SYSTEM=Linux-egcs
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4STTDIR=/afs/cern.ch/sw/geant4/stt/$REF/testtools/geant4/tests/tools
  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4LIB=$G4WORKDIR/lib
  export CLHEP_BASE_DIR=/afs/cern.ch/sw/geant4/dev/CLHEP/Linux-egcs/$CLHEP_VERSION
#
# The G4SYSTEm (and directory) was changed!
#  export CLHEP_LIB=CLHEP-egcs

  export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${CLHEP_BASE_DIR}/lib:${LD_LIBRARY_PATH}
  
  # G4 build flags :
  ######export G4UI_BUILD_XM_SESSION=1
#  export G4VIS_BUILD_OPENGLXM_DRIVER=1
#  export G4VIS_BUILD_OPENGLX_DRIVER=1
#  export G4VIS_BUILD_OIX_DRIVER=1
fi

if [ $UNAMEN = pcg4speed.cern.ch ]; then
  export DEBOPT=${DEBOPT}_slc3_icc91
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
  export G4SYSTEM=Linux-icc
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4STTDIR=/afs/cern.ch/sw/geant4/stt/$REF/testtools/geant4/tests/tools
  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4LIB=$G4WORKDIR/lib
#  export CLHEP_BASE_DIR=/afs/cern.ch/sw/geant4/dev/CLHEP/Linux-g++/$CLHEP_VERSION
  export CLHEP_BASE_DIR=$G4WORKDIR/clhep

  export XERCESCROOT=/afs/cern.ch/sw/geant4/dev/XercesC/2.8.0/slc3_icc91

  export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$XERCESCROOT/lib


# To link with AIDA (paths to libraries)
#export G4ANALYSIS_AIDA_CONFIG_LIBS="-L/afs/cern.ch/sw/lhcxx/specific/redhat61/gcc-2.95.2/5.0.1/lib -lAnaphe_AIDA_AnalysisFactory_native -lAnaphe_AIDA_Annotation_native -lAnaphe_AIDA_Histogram_native -lAnaphe_AIDA_Tree_native -lAnaphe_AIDA_HBookStore -lg2c-forMinuit -lHepUtilities -lCLHEP -lAnaphe_AIDA_Tuple_native"

#export G4ANALYSIS_AIDA_CONFIG_CFLAGS="-I/afs/cern.ch/sw/contrib/AIDA/3.0/src/cpp -I/afs/cern.ch/sw/lhcxx/specific/redhat61/gcc-2.95.2/5.0.1/include"

#export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:"/afs/cern.ch/sw/lhcxx/specific/redhat61/gcc-2.95.2/5.0.1/lib"

. /afs/cern.ch/sw/lcg/external/icc/9.1.039/slc3_ia32_gcc323/bin/iccvars.sh
export INTEL_LICENSE_FILE="1988@licman1.cern.ch"

# Only for one advanced example: 
export CCAL_CONFPATH=./dataconf
export CCAL_SENSITIVECONF=g4testbeamhcal96.conf
export CCAL_GEOMETRYCONF=testbeamhcal96.conf
export CCAL_GLOBALPATH=./dataglobal
export CCAL_GEOMPATH=./datageom
export CCAL_VISPATH=./datavis


  # Shareable library
  #####################
#  export G4LIB_BUILD_SHARED=1
#  export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${LD_LIBRARY_PATH}
  
  # G4 build flags :
  ######export G4UI_BUILD_XM_SESSION=1
#  export G4VIS_BUILD_OPENGLXM_DRIVER=1
#  export G4VIS_BUILD_OPENGLX_DRIVER=1
#  export G4VIS_BUILD_OIX_DRIVER=1
  export G4VIS_BUILD_OPENGLX_DRIVER=1
#  export G4VIS_USE_OPENGLX=1
fi


if [ $UNAMEN = pcg4speedBAC ]; then
  export DEBOPT=${DEBOPT}_NEWGCC
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
  export G4SYSTEM=Linux-g++
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4STTDIR=/afs/cern.ch/sw/geant4/stt/$REF/testtools/geant4/tests/tools
  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4LIB=$G4WORKDIR/lib
  export CLHEP_BASE_DIR=/afs/cern.ch/sw/geant4/dev/CLHEP/Linux-g++/$CLHEP_VERSION

# To link with AIDA (paths to libraries)
export G4ANALYSIS_AIDA_CONFIG_LIBS="-L/afs/cern.ch/sw/lhcxx/specific/redhat61/gcc-2.95.2/5.0.1/lib -lAnaphe_AIDA_AnalysisFactory_native -lAnaphe_AIDA_Annotation_native -lAnaphe_AIDA_Histogram_native -lAnaphe_AIDA_Tree_native -lAnaphe_AIDA_HBookStore -lg2c-forMinuit -lHepUtilities -lCLHEP -lAnaphe_AIDA_Tuple_native"

export G4ANALYSIS_AIDA_CONFIG_CFLAGS="-I/afs/cern.ch/sw/contrib/AIDA/3.0/src/cpp -I/afs/cern.ch/sw/lhcxx/specific/redhat61/gcc-2.95.2/5.0.1/include"

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:"/afs/cern.ch/sw/lhcxx/specific/redhat61/gcc-2.95.2/5.0.1/lib"

# Only for one advanced example: 
export CCAL_CONFPATH=./dataconf
export CCAL_SENSITIVECONF=g4testbeamhcal96.conf
export CCAL_GEOMETRYCONF=testbeamhcal96.conf
export CCAL_GLOBALPATH=./dataglobal
export CCAL_GEOMPATH=./datageom
export CCAL_VISPATH=./datavis


  # Shareable library
  #####################
#  export G4LIB_BUILD_SHARED=1
#  export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${LD_LIBRARY_PATH}
  
  # G4 build flags :
  ######export G4UI_BUILD_XM_SESSION=1
#  export G4VIS_BUILD_OPENGLXM_DRIVER=1
#  export G4VIS_BUILD_OPENGLX_DRIVER=1
#  export G4VIS_BUILD_OIX_DRIVER=1
fi


if [ $UNAMEN = pcitapi22 ]; then
  export DEBOPT=${DEBOPT}_NEWGCC
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
  export G4SYSTEM=Linux-g++
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4STTDIR=/afs/cern.ch/sw/geant4/stt/$REF/testtools/geant4/tests/tools
  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4LIB=$G4WORKDIR/lib
  export CLHEP_BASE_DIR=/afs/cern.ch/sw/geant4/dev/CLHEP/Linux-g++/$CLHEP_VERSION

  # Shareable library
  #####################
  export G4LIB_BUILD_SHARED=1
  export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${LD_LIBRARY_PATH}
  
  # G4 build flags :
  ######export G4UI_BUILD_XM_SESSION=1
#  export G4VIS_BUILD_OPENGLXM_DRIVER=1
#  export G4VIS_BUILD_OPENGLX_DRIVER=1
#  export G4VIS_BUILD_OIX_DRIVER=1
fi

if [ $UNAMEN = macphsft02.cern.ch ]; then
  export DEBOPT=${DEBOPT}_104_386_g401
  export CVSROOT=:ext:stesting@Geant4.cvs.cern.ch:/cvs/Geant4
  export CVS_RSH=ssh
  export G4SYSTEM=Darwin-g++
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4STTDIR=/afs/cern.ch/sw/geant4/stt/$REF/testtools/geant4/tests/tools
  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4LIB=$G4WORKDIR/lib
# Take CLHEP with links to lcg area
  export CLHEP_BASE_DIR=$G4WORKDIR/clhep
  
  export XERCESCROOT=/afs/cern.ch/sw/geant4/dev/XercesC/2.8.0/macosx104_gcc40
  
#  export G4_NO_CBRT=1

## Compiler
############
# Get gfortran
  export PATH=$PATH:/usr/local/bin
  
  # Shareable library
  #####################
#  export G4LIB_BUILD_SHARED=1
# export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${CLHEP_BASE_DIR}/lib:${LD_LIBRARY_PATH}
  export LD_LIBRARY_PATH=$XERCESCROOT/lib:/usr/local/lib
#  export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${LD_LIBRARY_PATH}

  # G4 build flags :
  ######export G4UI_BUILD_XM_SESSION=1
#  export G4VIS_BUILD_OPENGLXM_DRIVER=1
  export G4VIS_BUILD_OPENGLX_DRIVER=1
  export G4VIS_USE_OPENGLX=1
#  export G4VIS_BUILD_OIX_DRIVER=1
fi

if [ $UNAMEN = pcgeant2.cern.ch ]; then
  export DEBOPT=${DEBOPT}_slc3_323
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
  export G4SYSTEM=Linux-g++
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4STTDIR=/afs/cern.ch/sw/geant4/stt/$REF/testtools/geant4/tests/tools
  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/Linux-g++/$DEBOPT
  export G4LIB=$G4WORKDIR/lib
# Take CLHEP with links to lcg area
  export CLHEP_BASE_DIR=$G4WORKDIR/clhep 

  export XERCESCROOT=/afs/cern.ch/sw/geant4/dev/XercesC/2.8.0/slc3_gcc323

  export LD_LIBRARY_PATH=$XERCESCROOT/lib

#  export G4_NO_CBRT=1

## Compiler
############

  # Shareable library
  #####################
#  export G4LIB_BUILD_SHARED=1
# export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${CLHEP_BASE_DIR}/lib:${LD_LIBRARY_PATH}
#  export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${LD_LIBRARY_PATH}

  # G4 build flags :
  ######export G4UI_BUILD_XM_SESSION=1
#  export G4VIS_BUILD_OPENGLXM_DRIVER=1
  export G4VIS_BUILD_OPENGLX_DRIVER=1
  export G4VIS_USE_OPENGLX=1
#  export G4VIS_BUILD_OIX_DRIVER=1
fi

if [ $UNAMEN = lxbuild021.cern.ch ]; then
  export DEBOPT=${DEBOPT}_slc4_i32
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
  export G4SYSTEM=Linux-g++
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4STTDIR=/afs/cern.ch/sw/geant4/stt/$REF/testtools/geant4/tests/tools
  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/Linux-g++/$DEBOPT
  export G4LIB=$G4WORKDIR/lib

# Take CLHEP with links to lcg area
  export CLHEP_BASE_DIR=$G4WORKDIR/clhep
  
  export XERCESCROOT=/afs/cern.ch/sw/geant4/dev/XercesC/2.8.0/slc4_gcc346
  
#  export CLHEP_BASE_DIR=/scratch/stesting/clhep
  
  # Shareable library
  #####################
#  export G4LIB_BUILD_SHARED=1
#  export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${CLHEP_BASE_DIR}/lib:${LD_LIBRARY_PATH}
  export LD_LIBRARY_PATH=$XERCESCROOT/lib
#  export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${LD_LIBRARY_PATH}

  # G4 build flags :
  ######export G4UI_BUILD_XM_SESSION=1
#  export G4VIS_BUILD_OPENGLXM_DRIVER=1
  export G4VIS_BUILD_OPENGLX_DRIVER=1
  export G4VIS_USE_OPENGLX=1
#  export G4VIS_BUILD_OIX_DRIVER=1
fi

if [ $UNAMEN = pcg4speed2.cern.ch ]; then
  export DEBOPT=${DEBOPT}_slc4_i32
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
  export G4SYSTEM=Linux-g++
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4STTDIR=/afs/cern.ch/sw/geant4/stt/$REF/testtools/geant4/tests/tools
  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/Linux-g++/$DEBOPT
  export G4LIB=$G4WORKDIR/lib

# Take CLHEP with links to lcg area
  export CLHEP_BASE_DIR=$G4WORKDIR/clhep
#  export CLHEP_BASE_DIR=/scratch/stesting/clhep
  
  # Shareable library
  #####################
#  export G4LIB_BUILD_SHARED=1
  export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${CLHEP_BASE_DIR}/lib:${LD_LIBRARY_PATH}
#  export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${LD_LIBRARY_PATH}

  # G4 build flags :
  ######export G4UI_BUILD_XM_SESSION=1
#  export G4VIS_BUILD_OPENGLXM_DRIVER=1
  export G4VIS_BUILD_OPENGLX_DRIVER=1
  export G4VIS_USE_OPENGLX=1
#  export G4VIS_BUILD_OIX_DRIVER=1
fi

if [ $UNAMEN = lxplus071.cern.ch ]; then
# Special G4 flag for FPE exceptions on Linux systems
 export G4FPE_DEBUG=1
 
  export DEBOPT=${DEBOPT}_slc4_i32
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
  export G4SYSTEM=Linux-g++
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4STTDIR=/afs/cern.ch/sw/geant4/stt/$REF/testtools/geant4/tests/tools
  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/Linux-g++/$DEBOPT
  export G4LIB=$G4WORKDIR/lib

# Take CLHEP with links to lcg area
  export CLHEP_BASE_DIR=$G4WORKDIR/clhep
  
  export XERCESCROOT=/afs/cern.ch/sw/geant4/dev/XercesC/2.8.0/slc4_gcc346

  # Shareable library
  #####################
#  export G4LIB_BUILD_SHARED=1
#  export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${CLHEP_BASE_DIR}/lib:${LD_LIBRARY_PATH}
  export LD_LIBRARY_PATH=$XERCESCROOT/lib
#  export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${LD_LIBRARY_PATH}

  # G4 build flags :
  ######export G4UI_BUILD_XM_SESSION=1
#  export G4VIS_BUILD_OPENGLXM_DRIVER=1
  export G4VIS_BUILD_OPENGLX_DRIVER=1
  export G4VIS_USE_OPENGLX=1
#  export G4VIS_BUILD_OIX_DRIVER=1
fi

if [ $UNAMEN = lxplus072.cern.ch ]; then
# Special G4 flag for FPE exceptions on Linux systems
 export G4FPE_DEBUG=1
 
  export DEBOPT=${DEBOPT}_slc4_i32
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
  export G4SYSTEM=Linux-g++
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4STTDIR=/afs/cern.ch/sw/geant4/stt/$REF/testtools/geant4/tests/tools
  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/Linux-g++/$DEBOPT
  export G4LIB=$G4WORKDIR/lib

# Take CLHEP with links to lcg area
  export CLHEP_BASE_DIR=$G4WORKDIR/clhep
  
  export XERCESCROOT=/afs/cern.ch/sw/geant4/dev/XercesC/2.8.0/slc4_gcc346

  # Shareable library
  #####################
#  export G4LIB_BUILD_SHARED=1
#  export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${CLHEP_BASE_DIR}/lib:${LD_LIBRARY_PATH}
  export LD_LIBRARY_PATH=$XERCESCROOT/lib
#  export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${LD_LIBRARY_PATH}

  # G4 build flags :
  ######export G4UI_BUILD_XM_SESSION=1
#  export G4VIS_BUILD_OPENGLXM_DRIVER=1
  export G4VIS_BUILD_OPENGLX_DRIVER=1
  export G4VIS_USE_OPENGLX=1
#  export G4VIS_BUILD_OIX_DRIVER=1
fi


if [ $UNAMEN = lxcert-amd64 ]; then
  export DEBOPT=${DEBOPT}_slc4_amd64
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
  export G4SYSTEM=Linux-g++
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4STTDIR=/afs/cern.ch/sw/geant4/stt/$REF/testtools/geant4/tests/tools
  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/Linux-g++/$DEBOPT
  export G4LIB=$G4WORKDIR/lib

# Take CLHEP with links to lcg area
  export CLHEP_BASE_DIR=$G4WORKDIR/clhep
  
  export EXTRALIBS=" -L/usr/X11R6/lib64/ "

  # Shareable library
  #####################
#  export G4LIB_BUILD_SHARED=1
  export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${CLHEP_BASE_DIR}/lib:${LD_LIBRARY_PATH}
#  export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${LD_LIBRARY_PATH}

  # G4 build flags :
  ######export G4UI_BUILD_XM_SESSION=1
#  export G4VIS_BUILD_OPENGLXM_DRIVER=1
  export G4VIS_BUILD_OPENGLX_DRIVER=1
  export G4VIS_USE_OPENGLX=1
#  export G4VIS_BUILD_OIX_DRIVER=1
fi

if [ $UNAMEN = lxplus098.cern.ch ]; then
  export DEBOPT=${DEBOPT}_slc4_amd64
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
  export G4SYSTEM=Linux-g++
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4STTDIR=/afs/cern.ch/sw/geant4/stt/$REF/testtools/geant4/tests/tools
  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/Linux-g++/$DEBOPT
  export G4LIB=$G4WORKDIR/lib

# Take CLHEP with links to lcg area
  export CLHEP_BASE_DIR=$G4WORKDIR/clhep
  
  export EXTRALIBS=" -L/usr/X11R6/lib64/ "

  # Shareable library
  #####################
#  export G4LIB_BUILD_SHARED=1
  export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${CLHEP_BASE_DIR}/lib:${LD_LIBRARY_PATH}
#  export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${LD_LIBRARY_PATH}

  # G4 build flags :
  ######export G4UI_BUILD_XM_SESSION=1
#  export G4VIS_BUILD_OPENGLXM_DRIVER=1
  export G4VIS_BUILD_OPENGLX_DRIVER=1
  export G4VIS_USE_OPENGLX=1
#  export G4VIS_BUILD_OIX_DRIVER=1
fi

if [ $UNAMEN = lxplus099.cern.ch ]; then
  export DEBOPT=${DEBOPT}_slc4_amd64
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
  export G4SYSTEM=Linux-g++
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4STTDIR=/afs/cern.ch/sw/geant4/stt/$REF/testtools/geant4/tests/tools
  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/Linux-g++/$DEBOPT
  export G4LIB=$G4WORKDIR/lib

# Take CLHEP with links to lcg area
  export CLHEP_BASE_DIR=$G4WORKDIR/clhep
  
  export EXTRALIBS=" -L/usr/X11R6/lib64/ "

  # Shareable library
  #####################
#  export G4LIB_BUILD_SHARED=1
  export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${CLHEP_BASE_DIR}/lib:${LD_LIBRARY_PATH}
#  export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${LD_LIBRARY_PATH}

  # G4 build flags :
  ######export G4UI_BUILD_XM_SESSION=1
#  export G4VIS_BUILD_OPENGLXM_DRIVER=1
  export G4VIS_BUILD_OPENGLX_DRIVER=1
  export G4VIS_USE_OPENGLX=1
#  export G4VIS_BUILD_OIX_DRIVER=1
fi

if [ $UNAMEN = lxbuild056.cern.ch ]; then
  export DEBOPT=${DEBOPT}_slc4_amd64
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
  export G4SYSTEM=Linux-g++
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4STTDIR=/afs/cern.ch/sw/geant4/stt/$REF/testtools/geant4/tests/tools
  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/Linux-g++/$DEBOPT
  export G4LIB=$G4WORKDIR/lib

# Take CLHEP with links to lcg area
  export CLHEP_BASE_DIR=$G4WORKDIR/clhep
  
  export EXTRALIBS=" -L/usr/X11R6/lib64/ "

  # Shareable library
  #####################
#  export G4LIB_BUILD_SHARED=1
  export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${CLHEP_BASE_DIR}/lib:${LD_LIBRARY_PATH}
#  export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${LD_LIBRARY_PATH}

  # G4 build flags :
  ######export G4UI_BUILD_XM_SESSION=1
#  export G4VIS_BUILD_OPENGLXM_DRIVER=1
  export G4VIS_BUILD_OPENGLX_DRIVER=1
  export G4VIS_USE_OPENGLX=1
#  export G4VIS_BUILD_OIX_DRIVER=1
fi

if [ $UNAMEN = pcgeant3.cern.ch ]; then
  export DEBOPT=${DEBOPT}_slc4_amd64
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
  export G4SYSTEM=Linux-g++
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4STTDIR=/afs/cern.ch/sw/geant4/stt/$REF/testtools/geant4/tests/tools
  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/Linux-g++/$DEBOPT
  export G4LIB=$G4WORKDIR/lib

# Take CLHEP with links to lcg area
  export CLHEP_BASE_DIR=$G4WORKDIR/clhep
  
  export XERCESCROOT=/afs/cern.ch/sw/geant4/dev/XercesC/2.8.0/slc4_x86_64_gcc346
  
  export EXTRALIBS=" -L/usr/X11R6/lib64/ "

  # Shareable library
  #####################
#  export G4LIB_BUILD_SHARED=1
#  export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${CLHEP_BASE_DIR}/lib:${LD_LIBRARY_PATH}
  export LD_LIBRARY_PATH=$XERCESCROOT/lib
  
#  export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${LD_LIBRARY_PATH}

  # G4 build flags :
  ######export G4UI_BUILD_XM_SESSION=1
#  export G4VIS_BUILD_OPENGLXM_DRIVER=1
  export G4VIS_BUILD_OPENGLX_DRIVER=1
  export G4VIS_USE_OPENGLX=1
#  export G4VIS_BUILD_OIX_DRIVER=1
fi


if [ $UNAMEN = pcgeant5.cern.ch ]; then
  export DEBOPT=${DEBOPT}_7.3_3.2
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
#  export G4SYSTEM=Linux-gO2
  export G4SYSTEM=Linux-g++
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4STTDIR=/afs/cern.ch/sw/geant4/stt/$REF/testtools/geant4/tests/tools
#  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/Linux-g++/$DEBOPT
  export G4LIB=$G4WORKDIR/lib
# Take CLHEP with links to lcg area
  export CLHEP_BASE_DIR=$G4WORKDIR/clhep

#  export CLHEP_BASE_DIR=/afs/cern.ch/sw/geant4/dev/CLHEP/pro/Linux-g++
#  export CLHEP_BASE_DIR=/home2/stesting/newA_CLHEP

#  export CLHEP_BASE_DIR=/home2/stesting/local
#  export G4_NO_CBRT=1

## Compiler
############
#. /afs/cern.ch/sw/geant4/dev/scripts/gcc32.sh 
. /afs/cern.ch/sw/geant4/dev/scripts/gcc323.sh 

  # Shareable library
  #####################
#  export G4LIB_BUILD_SHARED=1
  export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${CLHEP_BASE_DIR}/lib:${LD_LIBRARY_PATH}
#  export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${LD_LIBRARY_PATH}

  # G4 build flags :
  ######export G4UI_BUILD_XM_SESSION=1
#  export G4VIS_BUILD_OPENGLXM_DRIVER=1
  export G4VIS_BUILD_OPENGLX_DRIVER=1
  export G4VIS_USE_OPENGLX=1
#  export G4VIS_BUILD_OIX_DRIVER=1
fi

if [ $UNAMEN = lxcert-i386.cern.ch ]; then
  export DEBOPT=${DEBOPT}_lxcert
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
#  export G4SYSTEM=Linux-gO2
  export G4SYSTEM=Linux-g++
echo $G4SYSTEM
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4STTDIR=/afs/cern.ch/sw/geant4/stt/$REF/testtools/geant4/tests/tools
#  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4WORKDIR=/scratch/stesting/work/geant4
  export G4LIB=$G4WORKDIR/lib
#  Take CLHEP with links to lcg area
  export CLHEP_BASE_DIR=$G4WORKDIR/clhep
  
#  export CLHEP_BASE_DIR=/afs/cern.ch/sw/geant4/dev/CLHEP/pro/Linux-g++

#  export G4_NO_CBRT=1

## Compiler
############
#. /afs/cern.ch/sw/geant4/dev/scripts/gcc32.sh 
#. /afs/cern.ch/sw/geant4/dev/scripts/gcc323.sh 

  # Shareable library
  #####################
#  export G4LIB_BUILD_SHARED=1
  export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${CLHEP_BASE_DIR}/lib:${LD_LIBRARY_PATH}
#  export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${LD_LIBRARY_PATH}

  # G4 build flags :
  ######export G4UI_BUILD_XM_SESSION=1
  export G4VIS_BUILD_OPENGLXM_DRIVER=1
  export G4VIS_BUILD_OPENGLX_DRIVER=1
#  export G4VIS_BUILD_OIX_DRIVER=1
  export G4VIS_BUILD_OPENGLX_DRIVER=1
#  export G4VIS_USE_OPENLX=1

fi


if [ $UNAMEN = pcgeant4.cern.ch ]; then
  export DEBOPT=${DEBOPT}_slc4_amd64
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
  export G4SYSTEM=Linux-g++
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4STTDIR=/afs/cern.ch/sw/geant4/stt/$REF/testtools/geant4/tests/tools
  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/Linux-g++/$DEBOPT
  export G4LIB=$G4WORKDIR/lib

# Take CLHEP with links to lcg area
  export CLHEP_BASE_DIR=$G4WORKDIR/clhep
  
  export XERCESCROOT=/afs/cern.ch/sw/geant4/dev/XercesC/2.8.0/slc4_x86_64_gcc346
  
  export EXTRALIBS=" -L/usr/X11R6/lib64/ "

  # Shareable library
  #####################
#  export G4LIB_BUILD_SHARED=1
  export LD_LIBRARY_PATH=$XERCESCROOT/lib
#  export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${LD_LIBRARY_PATH}

  # G4 build flags :
  ######export G4UI_BUILD_XM_SESSION=1
#  export G4VIS_BUILD_OPENGLXM_DRIVER=1
  export G4VIS_BUILD_OPENGLX_DRIVER=1
  export G4VIS_USE_OPENGLX=1
#  export G4VIS_BUILD_OIX_DRIVER=1
fi


# pcgeant4 "backup" for pcg4speed
if [ $UNAMEN = pcgeant4XXX ]; then
  export DEBOPT=${DEBOPT}_NEWGCC
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
  export G4SYSTEM=Linux-g++
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4STTDIR=/afs/cern.ch/sw/geant4/stt/$REF/testtools/geant4/tests/tools
  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4LIB=$G4WORKDIR/lib
  export CLHEP_BASE_DIR=/afs/cern.ch/sw/geant4/dev/CLHEP/Linux-g++/$CLHEP_VERSION

  # Shareable library
  #####################
  export G4LIB_BUILD_SHARED=1
  export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${LD_LIBRARY_PATH}
  
  # G4 build flags :
fi


if [ $UNAMEN = XXXlxplus073 ]; then
  export DEBOPT=${DEBOPT}_7.3.1_2.95.2
#  export DEBOPT=${DEBOPT}
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
  export G4SYSTEM=Linux-g++
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4STTDIR=/afs/cern.ch/sw/geant4/stt/$REF/testtools/geant4/tests/tools
  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4LIB=$G4WORKDIR/lib
#  export CLHEP_BASE_DIR=/home/stesting/CLHEP
#  export CLHEP_BASE_DIR=/afs/cern.ch/sw/lhcxx/specific/redhat61/gcc-2.95.2/CLHEP/1.7.0.0/
#  export CLHEP_BASE_DIR=/afs/cern.ch/sw/geant4/dev/CLHEP/Linux-g++/$CLHEP_VERSION
  export CLHEP_BASE_DIR=/afs/cern.ch/sw/geant4/dev/CLHEP/Linux-g++/RH72/1.8/
#  export G4_NO_CBRT=1

  # Shareable library
  #####################
  export G4LIB_BUILD_SHARED=1
  export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${CLHEP_BASE_DIR}/lib:${LD_LIBRARY_PATH}
#  export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${LD_LIBRARY_PATH}

  # G4 build flags :
  ######export G4UI_BUILD_XM_SESSION=1
#  export G4VIS_BUILD_OPENGLXM_DRIVER=1
#  export G4VIS_BUILD_OPENGLX_DRIVER=1
#  export G4VIS_BUILD_OIX_DRIVER=1

    # Set alternative g++ 2.95.2 compiler
    #####################################
export PATH=/usr/local/gcc-alt-2.95.2/bin:${PATH}
export LD_LIBRARY_PATH=/usr/local/gcc-alt-2.95.2/lib:${LD_LIBRARY_PATH}

fi


if [ $UNAMEN = lxbuild005 ]; then
  echo $UNAMEN only prepared to manipulate source code
  export DEBOPT=${DEBOPT}_7.3.1_2.95.2
#  export DEBOPT=${DEBOPT}
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
  export G4SYSTEM=Linux-g++
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4STTDIR=/afs/cern.ch/sw/geant4/stt/$REF/testtools/geant4/tests/tools
  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4LIB=$G4WORKDIR/lib
#  export CLHEP_BASE_DIR=/home/stesting/CLHEP
#  export CLHEP_BASE_DIR=/afs/cern.ch/sw/lhcxx/specific/redhat61/gcc-2.95.2/CLHEP/1.7.0.0/
#  export CLHEP_BASE_DIR=/afs/cern.ch/sw/geant4/dev/CLHEP/Linux-g++/$CLHEP_VERSION
  export CLHEP_BASE_DIR=/afs/cern.ch/sw/geant4/dev/CLHEP/Linux-g++/RH72/1.8/
#  export G4_NO_CBRT=1

  # Shareable library
  #####################
  export G4LIB_BUILD_SHARED=1
  export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${CLHEP_BASE_DIR}/lib:${LD_LIBRARY_PATH}
#  export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${LD_LIBRARY_PATH}

  # G4 build flags :
  ######export G4UI_BUILD_XM_SESSION=1
#  export G4VIS_BUILD_OPENGLXM_DRIVER=1
#  export G4VIS_BUILD_OPENGLX_DRIVER=1
#  export G4VIS_BUILD_OIX_DRIVER=1

    # Set alternative g++ 2.95.2 compiler
    #####################################
export PATH=/usr/local/gcc-alt-2.95.2/bin:${PATH}
export LD_LIBRARY_PATH=/usr/local/gcc-alt-2.95.2/lib:${LD_LIBRARY_PATH}

fi








if [ $UNAMEN = XXXlxplus075 ]; then
  export DEBOPT=${DEBOPT}_7.3.1_2.95.2
#  export DEBOPT=${DEBOPT}
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
  export G4SYSTEM=Linux-g++
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4STTDIR=/afs/cern.ch/sw/geant4/stt/$REF/testtools/geant4/tests/tools
  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4LIB=$G4WORKDIR/lib
#  export CLHEP_BASE_DIR=/home/stesting/CLHEP
#  export CLHEP_BASE_DIR=/afs/cern.ch/sw/lhcxx/specific/redhat61/gcc-2.95.2/CLHEP/1.7.0.0/
#  export CLHEP_BASE_DIR=/afs/cern.ch/sw/geant4/dev/CLHEP/Linux-g++/$CLHEP_VERSION
  export CLHEP_BASE_DIR=/afs/cern.ch/sw/geant4/dev/CLHEP/Linux-g++/RH72/1.8/
#  export G4_NO_CBRT=1

  # Shareable library
  #####################
  export G4LIB_BUILD_SHARED=1
  export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${CLHEP_BASE_DIR}/lib:${LD_LIBRARY_PATH}
#  export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${LD_LIBRARY_PATH}

  # G4 build flags :
  ######export G4UI_BUILD_XM_SESSION=1
#  export G4VIS_BUILD_OPENGLXM_DRIVER=1
#  export G4VIS_BUILD_OPENGLX_DRIVER=1
#  export G4VIS_BUILD_OIX_DRIVER=1

    # Set alternative g++ 2.95.2 compiler
    #####################################
export PATH=/usr/local/gcc-alt-2.95.2/bin:${PATH}
export LD_LIBRARY_PATH=/usr/local/gcc-alt-2.95.2/lib:${LD_LIBRARY_PATH}

fi

if [ $UNAMEN = oplapro01 ]; then
  export DEBOPT=${DEBOPT}_7.3.1_2.95.2
#  export DEBOPT=${DEBOPT}
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
  export G4SYSTEM=Linux-icc
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4STTDIR=/afs/cern.ch/sw/geant4/stt/$REF/testtools/geant4/tests/tools
#  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4WORKDIR=/data2/stesting/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4LIB=$G4WORKDIR/lib
#  export CLHEP_BASE_DIR=/home/stesting/CLHEP
#  export CLHEP_BASE_DIR=/afs/cern.ch/sw/lhcxx/specific/redhat61/gcc-2.95.2/CLHEP/1.7.0.0/
#  export CLHEP_BASE_DIR=/afs/cern.ch/sw/geant4/dev/CLHEP/Linux-g++/$CLHEP_VERSION
#  export CLHEP_BASE_DIR=/afs/cern.ch/sw/lhcxx/specific/ia64/gcc-3.2/CLHEP/1.8.0.0
  export CLHEP_BASE_DIR=/data2/stesting/CLHEP/new_icc/1.8
#  export G4_NO_CBRT=1

  # Shareable library
  #####################
#  export G4LIB_BUILD_SHARED=1
  export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${CLHEP_BASE_DIR}/lib:${LD_LIBRARY_PATH}
#  export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${LD_LIBRARY_PATH}

  # G4 build flags :
  ######export G4UI_BUILD_XM_SESSION=1
#  export G4VIS_BUILD_OPENGLXM_DRIVER=1
#  export G4VIS_BUILD_OPENGLX_DRIVER=1
#  export G4VIS_BUILD_OIX_DRIVER=1

    # Set alternative g++ 3.2 compiler
    #####################################
#. /afs/cern.ch/sw/lcg/external/icc/7.1.006/rh73/opt/intel/compiler70/ia32/bin/iccvars.sh


# new? /afs/cern.ch/sw/lcg/external/icc/8.0/rh73_ia64/bin/
. /afs/cern.ch/sw/lcg/external/icc/8.0.055/intel_cc_80/bin/iccvars.sh
#export PATH=/opt/gcc-3.2/bin:${PATH}
#export LD_LIBRARY_PATH=/opt/gcc-3.2/lib:${LD_LIBRARY_PATH}

fi

if [ $UNAMEN = OLDoplapro01 ]; then
  export DEBOPT=${DEBOPT}_7.3.1_2.95.2
#  export DEBOPT=${DEBOPT}
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
  export G4SYSTEM=Linux-g++
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4STTDIR=/afs/cern.ch/sw/geant4/stt/$REF/testtools/geant4/tests/tools
  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4LIB=$G4WORKDIR/lib
#  export CLHEP_BASE_DIR=/home/stesting/CLHEP
#  export CLHEP_BASE_DIR=/afs/cern.ch/sw/lhcxx/specific/redhat61/gcc-2.95.2/CLHEP/1.7.0.0/
#  export CLHEP_BASE_DIR=/afs/cern.ch/sw/geant4/dev/CLHEP/Linux-g++/$CLHEP_VERSION
  export CLHEP_BASE_DIR=/afs/cern.ch/sw/lhcxx/specific/ia64/gcc-3.2/CLHEP/1.8.0.0
#  export G4_NO_CBRT=1

  # Shareable library
  #####################
  export G4LIB_BUILD_SHARED=1
  export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${CLHEP_BASE_DIR}/lib:${LD_LIBRARY_PATH}
#  export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${LD_LIBRARY_PATH}

  # G4 build flags :
  ######export G4UI_BUILD_XM_SESSION=1
#  export G4VIS_BUILD_OPENGLXM_DRIVER=1
#  export G4VIS_BUILD_OPENGLX_DRIVER=1
#  export G4VIS_BUILD_OIX_DRIVER=1

    # Set alternative g++ 3.2 compiler
    #####################################
export PATH=/opt/gcc-3.2/bin:${PATH}
export LD_LIBRARY_PATH=/opt/gcc-3.2/lib:${LD_LIBRARY_PATH}

fi

if [ $UNAMEN = tbed0087 ]; then
  export DEBOPT=${DEBOPT}_AFS
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
  export G4SYSTEM=Linux-g++
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4STTDIR=/afs/cern.ch/sw/geant4/stt/$REF/testtools/geant4/tests/tools
  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4LIB=$G4WORKDIR/lib
#  export CLHEP_BASE_DIR=/home/stesting/CLHEP
#  export CLHEP_BASE_DIR=/afs/cern.ch/sw/lhcxx/specific/redhat61/gcc-2.95.2/CLHEP/1.7.0.0/
  export CLHEP_BASE_DIR=/afs/cern.ch/sw/geant4/dev/CLHEP/Linux-g++/$CLHEP_VERSION
#  export G4_NO_CBRT=1

# ALternative compiler
#################################
export PATH=/usr/local/gcc-alt-2.95.2/bin:${PATH}
export LD_LIBRARY_PATH=/usr/local/gcc-alt-2.95.2/lib:${LD_LIBRARY_PATH}

  # Shareable library
  #####################
  export G4LIB_BUILD_SHARED=1
  export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${CLHEP_BASE_DIR}/lib:${LD_LIBRARY_PATH}
#  export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${LD_LIBRARY_PATH}

  # G4 build flags :
  ######export G4UI_BUILD_XM_SESSION=1
#  export G4VIS_BUILD_OPENGLXM_DRIVER=1
#  export G4VIS_BUILD_OPENGLX_DRIVER=1
#  export G4VIS_BUILD_OIX_DRIVER=1
fi

if [ $UNAMEN = tbed0088 ]; then
  export DEBOPT=${DEBOPT}_AFS
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
  export G4SYSTEM=Linux-g++
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4STTDIR=/afs/cern.ch/sw/geant4/stt/$REF/testtools/geant4/tests/tools
  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4LIB=$G4WORKDIR/lib
#  export CLHEP_BASE_DIR=/home/stesting/CLHEP
#  export CLHEP_BASE_DIR=/afs/cern.ch/sw/lhcxx/specific/redhat61/gcc-2.95.2/CLHEP/1.7.0.0/
  export CLHEP_BASE_DIR=/afs/cern.ch/sw/geant4/dev/CLHEP/Linux-g++/$CLHEP_VERSION
#  export G4_NO_CBRT=1

# ALternative compiler
#################################
export PATH=/usr/local/gcc-alt-2.95.2/bin:${PATH}
export LD_LIBRARY_PATH=/usr/local/gcc-alt-2.95.2/lib:${LD_LIBRARY_PATH}

  # Shareable library
  #####################
  export G4LIB_BUILD_SHARED=1
  export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${CLHEP_BASE_DIR}/lib:${LD_LIBRARY_PATH}
#  export LD_LIBRARY_PATH=$G4LIB/$G4SYSTEM:${LD_LIBRARY_PATH}

  # G4 build flags :
  ######export G4UI_BUILD_XM_SESSION=1
#  export G4VIS_BUILD_OPENGLXM_DRIVER=1
#  export G4VIS_BUILD_OPENGLX_DRIVER=1
#  export G4VIS_BUILD_OIX_DRIVER=1
fi

if [ `uname -n | grep sgmedia` ]; then
  export G4USE_OSPACE=1
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
  export G4SYSTEM=SGI-CC
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4STTDIR=/afs/cern.ch/sw/geant4/stt/$REF/testtools/geant4/tests/tools
  export G4LIB=$G4WORKDIR/lib
  # G4 build flags :
  ######export G4UI_BUILD_XM_SESSION=1
#  export G4VIS_BUILD_OPENGLXM_DRIVER=1
#  export G4VIS_BUILD_OPENGLX_DRIVER=1
#  export G4VIS_BUILD_OIX_DRIVER=1
fi


UNAMEN=`uname -n `
if [ X"$UNAMEN" = X"pcl3eth4"  ]; then
UNAME=`uname `
if [ "$UNAME" = "AIX" ]; then
 export G4SYSTEM=AIX-xlC
else
 export G4SYSTEM=Linux-g++
fi	
export G4VIS_DEBUG=1
if [ "$UNAME" = "Linux" ]; then
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
  export G4SYSTEM=Linux-g++
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4LIB=$G4WORKDIR/lib
export G4VIS_BUILD_OPENGLX_DRIVER=1
export G4VIS_USE_OPENGLX=1
export G4VIS_BUILD_OPENGLXM_DRIVER=1
export G4VIS_USE_OPENGLXM=1
export CLHEP_BASE_DIR=/afs/cern.ch/sw/geant4/dev/CLHEP/Linux-g++/$CLHEP_VERSION
#/usr/local/CLHEP1.3/CLHEP
export CLHEP_LIB=CLHEP
#CLHEP-c++
export RWBASE=/usr/local
export RWINC=/usr/local/include
export OGLHOME=/usr/local
fi
export XM_INSTALLED=1
export XKEYSYMDB=/usr/lib/X11/XKeysymDB
export G4UI_USE_TERMINAL=1
export G4UI_USE_GAG=1
export G4UI_BUILD_XM_SESSION=1
fi
#
if [ `uname -n` = aleph ] ; then
export G4USE_OSPACE=1
export CVSROOT=:pserver:barrand@g4cvs.cern.ch:/afs/cern.ch/sw/geant4/cvs
export G4INSTALL=/geant4/stt/$REF/src/geant4
export G4SYSTEM=HP-aCC
export G4WORKDIR=/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
export G4DEBUG=1
# G4 build flags :
export G4UI_BUILD_XM_SESSION=1
export G4VIS_BUILD_OPENGLXM_DRIVER=1
export G4VIS_BUILD_OPENGLX_DRIVER=1
export G4VIS_BUILD_OIX_DRIVER=1
# G4 use flags :
export G4UI_USE_XM=1
export G4VIS_USE_OPENGLXM=1
export G4VIS_USE_OPENGLX=1
export G4VIS_USE_OIX=1
# Specific :
export CLHEP_BASE_DIR=/geant4/HP-UX
export OGLHOME=/geant4/HP-UX
export OIVFLAGS="-I/geant4/HP-UX/include/SoFree"
export OIVLIBS="-L/geant4/HP-UX/lib -lhepvisXt -lhepvis -lSoFreeXt -lSoFree"
export SOFREEUSER=/projects/SoFree/user/
fi
#
if [ `uname -n` = asc ] ; then
export G4USE_OSPACE=1
export CVSROOT=':pserver:barrand@g4cvs.cern.ch:/afs/cern.ch/sw/geant4/cvs'
export G4INSTALL=/geant4/stt/$REF/src/geant4
export G4SYSTEM=OSF1
export G4WORKDIR=/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
export G4DEBUG=1
# G4 build flags :
export G4UI_BUILD_XM_SESSION=1
export G4UI_BUILD_XAW_SESSION=1
export G4VIS_BUILD_OPENGLXM_DRIVER=1
export G4VIS_BUILD_OPENGLX_DRIVER=1
export G4VIS_BUILD_OIX_DRIVER=1
# G4 use flags :
export G4UI_USE_XM=1
export G4UI_USE_XAW=1
export G4VIS_USE_OPENGLXM=1
export G4VIS_USE_OPENGLX=1
export G4VIS_USE_OIX=1
# Specific :
export RWBASE=/geant4/OSF1
export CLHEP_BASE_DIR=/geant4/OSF1
export OGLHOME=/geant4/OSF1
export OIVHOME=/geant4/OpenInventor2.4.1
export OIVFLAGS="-I$OIVHOME/include -I/geant4/OSF1/include"
export OIVLIBS="-L/geant4/OSF1/lib -lhepvisXt -lhepvis  -L$OIVHOME/lib -lInventorXt -lInventor -limage"
export XENVIRONMENT=$OIVHOME/app-defaults/Inventor
fi

# Global environment
