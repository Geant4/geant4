#
# This file is sourced by setup.sh script.
# 
####################################################
####################################################
### Specific STT member machines configurations  ###
####################################################
####################################################
#

if [ `pwd | grep ref+` ]; then
  REF=ref+
else
  REF=ref
fi

if [ $DEBUG ]; then
  DEBOPT=debug
else
  DEBOPT=optim
fi

UNAMEN=`uname -n `
ANS=`uname -n | grep rsplus`
if [ X`uname -n | grep rsplus` != X  -o "$UNAMEN" = "shift51" ]; then
  export CVSROOT=/afs/cern.ch/rd44/cvs
  export G4SYSTEM=AIX-xlC
  export G4INSTALL=/afs/cern.ch/rd44/stt/$REF/src/geant4
  export G4WORKDIR=/afs/cern.ch/rd44/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4LIB=$G4WORKDIR/lib
  # G4 build flags :
#  export G4UI_BUILD_TERMINAL_SESSION=1
#  export G4UI_BUILD_GAG_SESSION=1
  #####export G4UI_BUILD_XM_SESSION=1
  #####export G4VIS_BUILD_OPENGLXM_DRIVER=1
#  export G4VIS_BUILD_OPENGLX_DRIVER=1
#  export OGLHOME=/afs/cern.ch/rd44/dev/Mesa/Mesa-1.2.8
#  export G4VIS_BUILD_RAYX_DRIVER=1
#  export G4VIS_BUILD_RAYX_DRIVER=1
  #####export G4VIS_BUILD_OIX_DRIVER=1
#  export G4VIS_BUILD_DAWN_DRIVER=1
#  export DAWN_BSD_UNIX_DOMAIN=1
#  export DAWN_HOME=/afs/cern.ch/rd44/dev/DAWN/AIX-AFS
#  export G4VIS_BUILD_DAWNFILE_DRIVER=1
#  export G4VIS_BUILD_VRML_DRIVER=1
#  export G4VIS_BUILD_VRMLFILE_DRIVER=1
fi

if [ `uname -n | grep sunasd1` ]; then
  export CVSROOT=/afs/cern.ch/rd44/cvs
  export G4SYSTEM=SUN-CC
  export G4INSTALL=/afs/cern.ch/rd44/stt/$REF/src/geant4
  export G4WORKDIR=/afs/cern.ch/rd44/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4LIB=$G4WORKDIR/lib
  # G4 build flags :
#  export G4UI_BUILD_TERMINAL_SESSION=1
#  export G4UI_BUILD_GAG_SESSION=1
  #######export G4UI_BUILD_XM_SESSION=1
  #######export G4VIS_BUILD_OPENGLXM_DRIVER=1
  #######export G4VIS_BUILD_OPENGLX_DRIVER=1
#  export G4VIS_BUILD_RAYX_DRIVER=1
  #######export G4VIS_BUILD_OIX_DRIVER=1
  #######export G4VIS_BUILD_DAWN_DRIVER=1
#  export G4VIS_BUILD_DAWNFILE_DRIVER=1
  #######export G4VIS_BUILD_VRML_DRIVER=1
#  export G4VIS_BUILD_VRMLFILE_DRIVER=1
fi

if [ `uname -n | grep suncmsb` ]; then
  export CVSROOT=/afs/cern.ch/rd44/cvs
  export G4SYSTEM=SUN-CC
  export G4INSTALL=/afs/cern.ch/rd44/stt/$REF/src/geant4
  export G4WORKDIR=/afs/cern.ch/rd44/stt/$REF/$G4SYSTEM/$DEBOPT
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
#  export G4VIS_BUILD_RAYX_DRIVER=1
  #######export G4VIS_BUILD_OIX_DRIVER=1
  #######export G4VIS_BUILD_DAWN_DRIVER=1
#  export G4VIS_BUILD_DAWNFILE_DRIVER=1
  #######export G4VIS_BUILD_VRML_DRIVER=1
#  export G4VIS_BUILD_VRMLFILE_DRIVER=1
fi

if [ `uname -n | grep sungeant` ]; then
  export CVSROOT=/afs/cern.ch/rd44/cvs
  export G4SYSTEM=SUN-CC
  export G4INSTALL=/afs/cern.ch/rd44/stt/$REF/src/geant4
  export G4WORKDIR=/afs/cern.ch/rd44/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4LIB=$G4WORKDIR/lib
  # G4 build flags :
#  export G4UI_BUILD_TERMINAL_SESSION=1
#  export G4UI_BUILD_GAG_SESSION=1
  #######export G4UI_BUILD_XM_SESSION=1
  #######export G4VIS_BUILD_OPENGLXM_DRIVER=1
  #######export G4VIS_BUILD_OPENGLX_DRIVER=1
#  export G4VIS_BUILD_RAYX_DRIVER=1
  #######export G4VIS_BUILD_OIX_DRIVER=1
  #######export G4VIS_BUILD_DAWN_DRIVER=1
#  export G4VIS_BUILD_DAWNFILE_DRIVER=1
  #######export G4VIS_BUILD_VRML_DRIVER=1
#  export G4VIS_BUILD_VRMLFILE_DRIVER=1
fi


if [ `uname -n | grep hpplus` ]; then
  export CVSROOT=/afs/cern.ch/rd44/cvs
  export G4SYSTEM=HP-aCC
  export G4INSTALL=/afs/cern.ch/rd44/stt/$REF/src/geant4
  export G4WORKDIR=/afs/cern.ch/rd44/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4LIB=$G4WORKDIR/lib
  # G4 build flags :
#  export G4UI_BUILD_TERMINAL_SESSION=1
#  export G4UI_BUILD_GAG_SESSION=1
  ######export G4UI_BUILD_XM_SESSION=1
#  export G4VIS_BUILD_OPENGLXM_DRIVER=1
#  export G4VIS_BUILD_OPENGLX_DRIVER=1
#  export G4VIS_BUILD_RAYX_DRIVER=1
#  export G4VIS_BUILD_OIX_DRIVER=1
#  export G4VIS_BUILD_DAWN_DRIVER=1
#  export G4VIS_BUILD_DAWNFILE_DRIVER=1
#  export G4VIS_BUILD_VRML_DRIVER=1
#  export G4VIS_BUILD_VRMLFILE_DRIVER=1
fi

if [ `uname -n | grep axcnsi` ]; then
  export CVSROOT=/afs/cern.ch/rd44/cvs
  export G4SYSTEM=DEC-cxx
  export G4INSTALL=/afs/cern.ch/rd44/stt/$REF/src/geant4
  export G4WORKDIR=/afs/cern.ch/rd44/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4LIB=$G4WORKDIR/lib
  # G4 build flags :
#  export G4UI_BUILD_TERMINAL_SESSION=1
#  export G4UI_BUILD_GAG_SESSION=1
  ######export G4UI_BUILD_XM_SESSION=1
#  export G4VIS_BUILD_OPENGLXM_DRIVER=1
#  export G4VIS_BUILD_OPENGLX_DRIVER=1
#  export G4VIS_BUILD_RAYX_DRIVER=1
#  export G4VIS_BUILD_OIX_DRIVER=1
#  export G4VIS_BUILD_DAWN_DRIVER=1
#  export G4VIS_BUILD_DAWNFILE_DRIVER=1
#  export G4VIS_BUILD_VRML_DRIVER=1
#  export G4VIS_BUILD_VRMLFILE_DRIVER=1
fi

if [ `uname -n | grep dxplus` ]; then
  export CVSROOT=/afs/cern.ch/rd44/cvs
  export G4SYSTEM=DEC-cxx
  export G4INSTALL=/afs/cern.ch/rd44/stt/$REF/src/geant4
  export G4WORKDIR=/afs/cern.ch/rd44/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4LIB=$G4WORKDIR/lib
  # G4 build flags :
#  export G4UI_BUILD_TERMINAL_SESSION=1
#  export G4UI_BUILD_GAG_SESSION=1
  ######export G4UI_BUILD_XM_SESSION=1
#  export G4VIS_BUILD_OPENGLXM_DRIVER=1
#  export G4VIS_BUILD_OPENGLX_DRIVER=1
#  export G4VIS_BUILD_RAYX_DRIVER=1
#  export G4VIS_BUILD_OIX_DRIVER=1
#  export G4VIS_BUILD_DAWN_DRIVER=1
#  export G4VIS_BUILD_DAWNFILE_DRIVER=1
#  export G4VIS_BUILD_VRML_DRIVER=1
#  export G4VIS_BUILD_VRMLFILE_DRIVER=1
fi

if [ `uname -n | grep pcitasd04` ]; then
  export CVSROOT=/afs/cern.ch/rd44/cvs
  export G4SYSTEM=Linux-g++
  export G4INSTALL=/afs/cern.ch/rd44/stt/$REF/src/geant4
  export G4WORKDIR=/afs/cern.ch/rd44/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4LIB=$G4WORKDIR/lib
  # G4 build flags :
#  export G4UI_BUILD_TERMINAL_SESSION=1
#  export G4UI_BUILD_GAG_SESSION=1
  ######export G4UI_BUILD_XM_SESSION=1
#  export G4VIS_BUILD_OPENGLXM_DRIVER=1
#  export G4VIS_BUILD_OPENGLX_DRIVER=1
#  export G4VIS_BUILD_RAYX_DRIVER=1
#  export G4VIS_BUILD_OIX_DRIVER=1
#  export G4VIS_BUILD_DAWN_DRIVER=1
#  export G4VIS_BUILD_DAWNFILE_DRIVER=1
#  export G4VIS_BUILD_VRML_DRIVER=1
#  export G4VIS_BUILD_VRMLFILE_DRIVER=1
fi


if [ `uname -n | grep sgmedia` ]; then
  export CVSROOT=/afs/cern.ch/rd44/cvs
  export G4SYSTEM=SGI-CC
  export G4INSTALL=/afs/cern.ch/rd44/stt/$REF/src/geant4
  export G4WORKDIR=/afs/cern.ch/rd44/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4LIB=$G4WORKDIR/lib
  # G4 build flags :
#  export G4UI_BUILD_TERMINAL_SESSION=1
#  export G4UI_BUILD_GAG_SESSION=1
  ######export G4UI_BUILD_XM_SESSION=1
#  export G4VIS_BUILD_OPENGLXM_DRIVER=1
#  export G4VIS_BUILD_OPENGLX_DRIVER=1
#  export G4VIS_BUILD_RAYX_DRIVER=1
#  export G4VIS_BUILD_OIX_DRIVER=1
#  export G4VIS_BUILD_DAWN_DRIVER=1
#  export G4VIS_BUILD_DAWNFILE_DRIVER=1
#  export G4VIS_BUILD_VRML_DRIVER=1
#  export G4VIS_BUILD_VRMLFILE_DRIVER=1
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
  export CVSROOT=/afs/cern.ch/rd44/cvs
  export G4SYSTEM=Linux-g++
  export G4INSTALL=/afs/cern.ch/rd44/stt/$REF/src/geant4
  export G4WORKDIR=/afs/cern.ch/rd44/stt/$REF/$G4SYSTEM/$DEBOPT
  export G4LIB=$G4WORKDIR/lib
export G4VIS_BUILD_OPENGLX_DRIVER=1
export G4VIS_USE_OPENGLX=1
export G4VIS_BUILD_OPENGLXM_DRIVER=1
export G4VIS_USE_OPENGLXM=1
export CLHEP_BASE_DIR=/afs/cern.ch/rd44/stt/$REF/Linux-g++/CLHEP
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
export CVSROOT=:pserver:barrand@g4cvs.cern.ch:/afs/cern.ch/rd44/cvs
export G4INSTALL=/geant4/stt/$REF/src/geant4
export G4SYSTEM=HP-aCC
export G4WORKDIR=/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
export G4DEBUG=1
# G4 build flags :
export G4UI_BUILD_XM_SESSION=1
export G4VIS_BUILD_OPENGLXM_DRIVER=1
export G4VIS_BUILD_OPENGLX_DRIVER=1
export G4VIS_BUILD_RAYX_DRIVER=1
export G4VIS_BUILD_OIX_DRIVER=1
export G4VIS_BUILD_DAWN_DRIVER=1
export G4VIS_BUILD_DAWNFILE_DRIVER=1
export G4VIS_BUILD_VRML_DRIVER=1
export G4VIS_BUILD_VRMLFILE_DRIVER=1
# G4 use flags :
export G4UI_USE_XM=1
export G4VIS_USE_OPENGLXM=1
export G4VIS_USE_OPENGLX=1
export G4VIS_USE_RAYX=1
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
export CVSROOT=':pserver:barrand@g4cvs.cern.ch:/afs/cern.ch/rd44/cvs'
export G4INSTALL=/geant4/stt/$REF/src/geant4
export G4SYSTEM=OSF1
export G4WORKDIR=/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
export G4DEBUG=1
# G4 build flags :
export G4UI_BUILD_XM_SESSION=1
export G4UI_BUILD_XAW_SESSION=1
export G4VIS_BUILD_OPENGLXM_DRIVER=1
export G4VIS_BUILD_OPENGLX_DRIVER=1
export G4VIS_BUILD_RAYX_DRIVER=1
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
export G4VIS_USE_RAYX=1
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

UNAMEN=`uname -n`
if [ `echo $UNAMEN | grep wagtail` ] ; then
export G4INSTALL=~/GEANT4/stt/$REF/src/geant4
export G4SYSTEM=Linux-g++
export G4WORKDIR=/scratch/stt/$REF/$G4SYSTEM/$DEBOPT
export G4LIB=$G4WORKDIR/lib
# G4 build flags :
export G4VIS_BUILD_OPENGLX_DRIVER=1
export G4VIS_BUILD_RAYX_DRIVER=1
export G4VIS_BUILD_DAWN_DRIVER=1
export G4VIS_BUILD_DAWNFILE_DRIVER=1
export G4VIS_BUILD_VRML_DRIVER=1
export G4VIS_BUILD_VRMLFILE_DRIVER=1
# G4 use flags :
export G4VIS_USE_OPENGLX=1
export G4VIS_USE_RAYX=1
export G4VIS_USE_DAWN=1
export G4VIS_USE_DAWNFILE=1
export G4VIS_USE_VRML=1
export G4VIS_USE_VRMLFILE=1
# Specific :
export MOMOPATH=$G4INSTALL/environments/Momo/tcltk
export CVSROOT=:pserver:allison@g4cvs.cern.ch:/afs/cern.ch/rd44/cvs
export CVSEDITOR=emacs
export CLHEP_BASE_DIR=/usr/local
export RWBASE=/usr/local
export XKEYSYMDB=/usr/lib/X11/XKeysymDB
export OGLHOME=/usr/X11
fi

