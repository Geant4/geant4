#
# This file is sourced by setup.sh script.
# 
####################################################
####################################################
### Specific STT member machines configurations  ###
####################################################
####################################################
#

REF=undefined
if [ `pwd | grep /stt/dev1/` ]; then
  REF=dev1
fi
if [ `pwd | grep /stt/dev2/` ]; then
  REF=dev2
fi
if [ `pwd | grep /stt/prod/` ]; then
  REF=prod
fi


if [ $G4DEBUG ]; then
  DEBOPT=debug
else
  DEBOPT=optim
fi

UNAMEN=`uname -n `
ANS=`uname -n | grep rsplus`
if [ X`uname -n | grep rsplus` != X  -o "$UNAMEN" = "shift51" ]; then
  export G4USE_OSPACE=1
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
  export G4SYSTEM=AIX-xlC
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4LIB=$G4WORKDIR/lib
  # G4 build flags :
#  export G4UI_BUILD_TERMINAL_SESSION=1
#  export G4UI_BUILD_GAG_SESSION=1
  #####export G4UI_BUILD_XM_SESSION=1
  #####export G4VIS_BUILD_OPENGLXM_DRIVER=1
  export G4VIS_BUILD_OPENGLX_DRIVER=1
  export OGLHOME=/afs/cern.ch/sw/geant4/dev/Mesa/Mesa-2.5
  #####export G4VIS_BUILD_OIX_DRIVER=1
  export G4VIS_BUILD_DAWN_DRIVER=1
#  export DAWN_BSD_UNIX_DOMAIN=1
#  export DAWN_HOME=/afs/cern.ch/sw/geant4/dev/DAWN/AIX-AFS
  export G4VIS_BUILD_DAWNFILE_DRIVER=1
  export G4VIS_BUILD_VRML_DRIVER=1
  export G4VIS_BUILD_VRMLFILE_DRIVER=1
fi

if [ `uname -n | grep sunasd1` ]; then
  export G4USE_OSPACE=1
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
  export G4SYSTEM=SUN-CC
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4LIB=$G4WORKDIR/lib
  # G4 build flags :
#  export G4UI_BUILD_TERMINAL_SESSION=1
#  export G4UI_BUILD_GAG_SESSION=1
  #######export G4UI_BUILD_XM_SESSION=1
  #######export G4VIS_BUILD_OPENGLXM_DRIVER=1
  #######export G4VIS_BUILD_OPENGLX_DRIVER=1
  #######export G4VIS_BUILD_OIX_DRIVER=1
  export G4VIS_BUILD_DAWN_DRIVER=1
  export G4VIS_BUILD_DAWNFILE_DRIVER=1
  export G4VIS_BUILD_VRML_DRIVER=1
  export G4VIS_BUILD_VRMLFILE_DRIVER=1
fi

if [ `uname -n | grep suncmsb` ]; then
  export G4USE_OSPACE=1
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
  export G4SYSTEM=SUN-CC
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4LIB=$G4WORKDIR/lib
#
  export CLHEP_BASE_DIR=/afs/cern.ch/user/s/stesting/work/clhep
  export CLHEP_LIB=CLHEP-CC
  export RWBASE=/afs/cern.ch/user/s/stesting/work/rogue
  # G4 build flags :
#  export G4UI_BUILD_TERMINAL_SESSION=1
#  export G4UI_BUILD_GAG_SESSION=1
  #######export G4UI_BUILD_XM_SESSION=1
  #######export G4VIS_BUILD_OPENGLXM_DRIVER=1
  #######export G4VIS_BUILD_OPENGLX_DRIVER=1
  #######export G4VIS_BUILD_OIX_DRIVER=1
  export G4VIS_BUILD_DAWN_DRIVER=1
  export G4VIS_BUILD_DAWNFILE_DRIVER=1
  export G4VIS_BUILD_VRML_DRIVER=1
  export G4VIS_BUILD_VRMLFILE_DRIVER=1
fi

if [ `uname -n | grep sungeant` ]; then
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  if [ $G4STTNONISO ]; then
    export G4SYSTEM=SUN-CC
    export DEBOPT=${DEBOPT}_NONISO
    export G4USE_OSPACE=1
    export PATH=`echo $PATH | sed s/SUNWspro50/SUNWspro/`
    # Persistency...
    if [ X$G4USE_HEPODBMS = X ]; then
      . $G4INSTALL/examples/extended/persistency/PersistentEx01/g4odbms_setup.sh
      export G4EXAMPLE_FDID=207
    fi
  else
    export G4SYSTEM=SUN-CC5
    export DEBOPT=${DEBOPT}_ISO
    unset G4USE_OSPACE
    export PATH=`echo $PATH | sed s/SUNWspro50/SUNWspro/`
    export PATH=`echo $PATH | sed s/SUNWspro/SUNWspro50/`
    # No Persistency...
    unset G4USE_HEPODBMS
  fi
  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4LIB=$G4WORKDIR/lib
  # G4 build flags :
#  export G4UI_BUILD_TERMINAL_SESSION=1
#  export G4UI_BUILD_GAG_SESSION=1
  #######export G4UI_BUILD_XM_SESSION=1
  #######export G4VIS_BUILD_OPENGLXM_DRIVER=1
  export G4VIS_BUILD_OPENGLX_DRIVER=1
  export OGLHOME=/usr/local
  export OGLFLAGS="-I$OGLHOME/include"
  export OGLLIBS="-L$OGLHOME/lib -lMesaGLU -lMesaGL"
  #######export G4VIS_BUILD_OIX_DRIVER=1
  export G4VIS_BUILD_DAWN_DRIVER=1
  export G4VIS_BUILD_DAWNFILE_DRIVER=1
  export G4VIS_BUILD_VRML_DRIVER=1
  export G4VIS_BUILD_VRMLFILE_DRIVER=1
fi


if [ `uname -n | grep hp` ]; then
  export G4USE_OSPACE=1
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
  export G4SYSTEM=HP-aCC
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4LIB=$G4WORKDIR/lib
  # G4 build flags :
#  export G4UI_BUILD_TERMINAL_SESSION=1
#  export G4UI_BUILD_GAG_SESSION=1
  ######export G4UI_BUILD_XM_SESSION=1
#  export G4VIS_BUILD_OPENGLXM_DRIVER=1
#  export G4VIS_BUILD_OPENGLX_DRIVER=1
#  export G4VIS_BUILD_OIX_DRIVER=1
  export G4VIS_BUILD_DAWN_DRIVER=1
  export G4VIS_BUILD_DAWNFILE_DRIVER=1
  export G4VIS_BUILD_VRML_DRIVER=1
  export G4VIS_BUILD_VRMLFILE_DRIVER=1
fi

if [ `uname -n | grep axcnsi` ]; then
  export G4USE_OSPACE=1
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
  export G4SYSTEM=DEC-cxx
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4LIB=$G4WORKDIR/lib
  # G4 build flags :
#  export G4UI_BUILD_TERMINAL_SESSION=1
#  export G4UI_BUILD_GAG_SESSION=1
  ######export G4UI_BUILD_XM_SESSION=1
#  export G4VIS_BUILD_OPENGLXM_DRIVER=1
#  export G4VIS_BUILD_OPENGLX_DRIVER=1
#  export G4VIS_BUILD_OIX_DRIVER=1
  export G4VIS_BUILD_DAWN_DRIVER=1
  export G4VIS_BUILD_DAWNFILE_DRIVER=1
  export G4VIS_BUILD_VRML_DRIVER=1
  export G4VIS_BUILD_VRMLFILE_DRIVER=1
fi

if [ X`uname -n | grep dxplus` != X  -o "$UNAMEN" = "dcosf01" ]; then
  if [ $G4STTNONISO ]; then
    export DEBOPT=${DEBOPT}_NONISO
    export G4USE_OSPACE=1
    export CLHEP_BASE_DIR=/afs/cern.ch/sw/geant4/dev/CLHEP/DEC-cxx/pro
  else
    export DEBOPT=${DEBOPT}_ISO
    unset G4USE_OSPACE
    export CLHEP_BASE_DIR=/afs/cern.ch/sw/geant4/dev/CLHEP/DEC-cxx/iso
  fi
  if [ "$UNAMEN" = "dcosf01" ]; then
    export CLHEP_LIB=CLHEP-cxx62
  fi
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
  export G4SYSTEM=DEC-cxx
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4LIB=$G4WORKDIR/lib
  # G4 build flags :
#  export G4UI_BUILD_TERMINAL_SESSION=1
#  export G4UI_BUILD_GAG_SESSION=1
  ######export G4UI_BUILD_XM_SESSION=1
#  export G4VIS_BUILD_OPENGLXM_DRIVER=1
#  export G4VIS_BUILD_OPENGLX_DRIVER=1
#  export G4VIS_BUILD_OIX_DRIVER=1
  export G4VIS_BUILD_DAWN_DRIVER=1
  export G4VIS_BUILD_DAWNFILE_DRIVER=1
  export G4VIS_BUILD_VRML_DRIVER=1
  export G4VIS_BUILD_VRMLFILE_DRIVER=1
fi

if [ `uname -n | grep pcitasd04` ]; then
  export G4USE_STL=1
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
  export G4SYSTEM=Linux-g++
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4LIB=$G4WORKDIR/lib
  # G4 build flags :
#  export G4UI_BUILD_TERMINAL_SESSION=1
#  export G4UI_BUILD_GAG_SESSION=1
  ######export G4UI_BUILD_XM_SESSION=1
#  export G4VIS_BUILD_OPENGLXM_DRIVER=1
#  export G4VIS_BUILD_OPENGLX_DRIVER=1
#  export G4VIS_BUILD_OIX_DRIVER=1
  export G4VIS_BUILD_DAWN_DRIVER=1
  export G4VIS_BUILD_DAWNFILE_DRIVER=1
  export G4VIS_BUILD_VRML_DRIVER=1
  export G4VIS_BUILD_VRMLFILE_DRIVER=1
fi

UNAMEN=`uname -n `
echo UNAMEN $UNAMEN
if [ $UNAMEN = pcgeant -o $UNAMEN = pcg4speed ]; then
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
  export G4SYSTEM=Linux-g++
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4LIB=$G4WORKDIR/lib
  # G4 build flags :
#  export G4UI_BUILD_TERMINAL_SESSION=1
#  export G4UI_BUILD_GAG_SESSION=1
  ######export G4UI_BUILD_XM_SESSION=1
#  export G4VIS_BUILD_OPENGLXM_DRIVER=1
#  export G4VIS_BUILD_OPENGLX_DRIVER=1
#  export G4VIS_BUILD_OIX_DRIVER=1
  export G4VIS_BUILD_DAWN_DRIVER=1
  export G4VIS_BUILD_DAWNFILE_DRIVER=1
  export G4VIS_BUILD_VRML_DRIVER=1
  export G4VIS_BUILD_VRMLFILE_DRIVER=1
fi


if [ `uname -n | grep sgmedia` ]; then
  export G4USE_OSPACE=1
  export CVSROOT=/afs/cern.ch/sw/geant4/cvs
  export G4SYSTEM=SGI-CC
  export G4INSTALL=/afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  export G4WORKDIR=/afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4LIB=$G4WORKDIR/lib
  # G4 build flags :
#  export G4UI_BUILD_TERMINAL_SESSION=1
#  export G4UI_BUILD_GAG_SESSION=1
  ######export G4UI_BUILD_XM_SESSION=1
#  export G4VIS_BUILD_OPENGLXM_DRIVER=1
#  export G4VIS_BUILD_OPENGLX_DRIVER=1
#  export G4VIS_BUILD_OIX_DRIVER=1
  export G4VIS_BUILD_DAWN_DRIVER=1
  export G4VIS_BUILD_DAWNFILE_DRIVER=1
  export G4VIS_BUILD_VRML_DRIVER=1
  export G4VIS_BUILD_VRMLFILE_DRIVER=1
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
export G4VIS_BUILD_DAWN_DRIVER=1
export DAWN_BSD_UNIX_DOMAIN=1
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
export CLHEP_BASE_DIR=/afs/cern.ch/sw/geant4/stt/$REF/Linux-g++/CLHEP
#/usr/local/CLHEP1.3/CLHEP
export CLHEP_LIB=CLHEP
#CLHEP-c++
export RWBASE=/usr/local
export RWINC=/usr/local/include
export OGLHOME=/usr/local
fi
export XM_INSTALLED=1
export XKEYSYMDB=/usr/lib/X11/XKeysymDB
export G4UI_BUILD_TERMINAL_SESSION=1
export G4UI_USE_TERMINAL=1
export G4UI_BUILD_GAG_SESSION=1
export G4UI_USE_GAG=1
export G4UI_BUILD_XM_SESSION=1
export G4VIS_BUILD_VRML_DRIVER=1
export G4VIS_BUILD_VRMLFILE_DRIVER=1
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
export G4VIS_BUILD_DAWN_DRIVER=1
export G4VIS_BUILD_DAWNFILE_DRIVER=1
export G4VIS_BUILD_VRML_DRIVER=1
export G4VIS_BUILD_VRMLFILE_DRIVER=1
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
export G4VIS_BUILD_DAWN_DRIVER=1
export G4VIS_BUILD_DAWNFILE_DRIVER=1
export G4VIS_BUILD_VRML_DRIVER=1
export G4VIS_BUILD_VRMLFILE_DRIVER=1
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
