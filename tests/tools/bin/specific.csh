
#
# This file is sourced by setup.csh script.
# 
####################################################
####################################################
### Specific STT member machines configurations  ###
####################################################
####################################################
#
####################################################
# Guy Barrand barrand@lal.in2p3.fr
####################################################

if ( `pwd | grep ref+` != "" ) then
  setenv REF ref+
else
  setenv REF ref
endif

if ( $G4DEBUG != "" ) then
  setenv DEBOPT debug
else
  setenv DEBOPT optim
endif

if ( `uname -n | grep rsplus` != "" ) then
  setenv CVSROOT /afs/cern.ch/sw/geant4/cvs
  setenv G4SYSTEM AIX-xlC
  setenv G4INSTALL /afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  setenv G4WORKDIR  /afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  setenv G4LIB $G4WORKDIR/lib
  # G4 build flags :
  setenv G4UI_BUILD_TERMINAL_SESSION 1
  setenv G4UI_BUILD_GAG_SESSION      1
  #####setenv G4UI_BUILD_XM_SESSION       1
  #####setenv G4VIS_BUILD_OPENGLXM_DRIVER 1
  setenv G4VIS_BUILD_OPENGLX_DRIVER  1
  setenv XKEYSYMDB /usr/lib/X11/XKeysymDB
  setenv OGLHOME /afs/cern.ch/sw/geant4/dev/Mesa/Mesa-1.2.8
  setenv G4VIS_BUILD_RAYX_DRIVER     1
  ##### setenv G4VIS_BUILD_OIX_DRIVER      1
  setenv G4VIS_BUILD_DAWN_DRIVER     1
  setenv DAWN_BSD_UNIX_DOMAIN 1
  setenv DAWN_HOME /afs/cern.ch/sw/geant4/dev/DAWN/AIX-AFS
  setenv G4VIS_BUILD_DAWNFILE_DRIVER 1
  setenv G4VIS_BUILD_VRML_DRIVER     1
  setenv G4VIS_BUILD_VRMLFILE_DRIVER 1
endif

if ( `uname -n | grep dxplus` != "" ) then
  setenv CVSROOT /afs/cern.ch/sw/geant4/cvs
  setenv G4SYSTEM DEC-cxx
  setenv G4INSTALL /afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  setenv G4WORKDIR  /afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  setenv G4LIB $G4WORKDIR/lib
  # G4 build flags :
  setenv G4UI_BUILD_TERMINAL_SESSION 1
  setenv G4UI_BUILD_GAG_SESSION      1
  #####setenv G4UI_BUILD_XM_SESSION       1
  #####setenv G4VIS_BUILD_OPENGLXM_DRIVER 1
  setenv G4VIS_BUILD_OPENGLX_DRIVER  1
  setenv XKEYSYMDB /usr/lib/X11/XKeysymDB
  setenv OGLHOME /afs/cern.ch/sw/geant4/dev/Mesa/Mesa-1.2.8
  setenv G4VIS_BUILD_RAYX_DRIVER     1
  ##### setenv G4VIS_BUILD_OIX_DRIVER      1
  setenv G4VIS_BUILD_DAWN_DRIVER     1
  setenv DAWN_BSD_UNIX_DOMAIN 1
  setenv DAWN_HOME /afs/cern.ch/sw/geant4/dev/DAWN/AIX-AFS
  setenv G4VIS_BUILD_DAWNFILE_DRIVER 1
  setenv G4VIS_BUILD_VRML_DRIVER     1
  setenv G4VIS_BUILD_VRMLFILE_DRIVER 1
endif

if ( `uname -n | grep sgmedia` != "" ) then
  setenv CVSROOT /afs/cern.ch/sw/geant4/cvs
  setenv G4SYSTEM SGI-CC
  setenv G4INSTALL /afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  setenv G4WORKDIR  /afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  setenv G4LIB $G4WORKDIR/lib
  # G4 build flags :
  #######setenv G4UI_BUILD_TERMINAL_SESSION 1
  #######setenv G4UI_BUILD_GAG_SESSION      1
  #######setenv G4UI_BUILD_XM_SESSION       1
  #######setenv G4VIS_BUILD_OPENGLXM_DRIVER 1
  #######setenv G4VIS_BUILD_OPENGLX_DRIVER  1
  #######setenv G4VIS_BUILD_RAYX_DRIVER     1
  #######setenv G4VIS_BUILD_OIX_DRIVER      1
  #######setenv G4VIS_BUILD_DAWN_DRIVER     1
  #######setenv G4VIS_BUILD_DAWNFILE_DRIVER 1
  #######setenv G4VIS_BUILD_VRML_DRIVER     1
  #######setenv G4VIS_BUILD_VRMLFILE_DRIVER 1
endif

if ( `uname -n | grep sun` != "" ) then
  setenv CVSROOT /afs/cern.ch/sw/geant4/cvs
  setenv G4SYSTEM SUN-CC
  setenv G4INSTALL /afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  setenv G4WORKDIR  /afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  setenv G4LIB $G4WORKDIR/lib
  # G4 build flags :
  setenv G4UI_BUILD_TERMINAL_SESSION 1
  setenv G4UI_BUILD_GAG_SESSION      1
  #######setenv G4UI_BUILD_XM_SESSION       1
  #######setenv G4VIS_BUILD_OPENGLXM_DRIVER 1
  #######setenv G4VIS_BUILD_OPENGLX_DRIVER  1
  setenv G4VIS_BUILD_RAYX_DRIVER     1
  #######setenv G4VIS_BUILD_OIX_DRIVER      1
  #######setenv G4VIS_BUILD_DAWN_DRIVER     1
  setenv G4VIS_BUILD_DAWNFILE_DRIVER 1
  #######setenv G4VIS_BUILD_VRML_DRIVER     1
  setenv G4VIS_BUILD_VRMLFILE_DRIVER 1
endif

if ( `uname -n | grep hpplus` != "" ) then
  setenv CVSROOT /afs/cern.ch/sw/geant4/cvs
  setenv G4SYSTEM HP-aCC
  setenv G4INSTALL /afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  setenv G4WORKDIR  /afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  setenv G4LIB $G4WORKDIR/lib
  # G4 build flags :
  setenv G4UI_BUILD_TERMINAL_SESSION 1
  setenv G4UI_BUILD_GAG_SESSION      1
  ######setenv G4UI_BUILD_XM_SESSION       1
  setenv G4VIS_BUILD_OPENGLXM_DRIVER 1
  setenv G4VIS_BUILD_OPENGLX_DRIVER  1
  setenv G4VIS_BUILD_RAYX_DRIVER     1
  setenv G4VIS_BUILD_OIX_DRIVER      1
  setenv G4VIS_BUILD_DAWN_DRIVER     1
  setenv G4VIS_BUILD_DAWNFILE_DRIVER 1
  setenv G4VIS_BUILD_VRML_DRIVER     1
  setenv G4VIS_BUILD_VRMLFILE_DRIVER 1
endif

if ( `uname -n` == aleph ) then
setenv CVSROOT :pserver:barrand@g4cvs.cern.ch:/afs/cern.ch/sw/geant4/cvs
setenv G4INSTALL                   /geant4/stt/$REF/src/geant4
setenv G4SYSTEM                    HP-aCC
setenv G4WORKDIR                   /geant4/stt/$REF/$G4SYSTEM/$DEBOPT
setenv G4DEBUG                     1
#setenv G4MAKESHLIB                 $G4INSTALL/config/makeshlib.sh
# G4 build flags :
setenv G4UI_BUILD_XM_SESSION       1
setenv G4VIS_BUILD_OPENGLXM_DRIVER 1
setenv G4VIS_BUILD_OPENGLX_DRIVER  1
setenv G4VIS_BUILD_RAYX_DRIVER     1
setenv G4VIS_BUILD_OIX_DRIVER      1
setenv G4VIS_BUILD_DAWN_DRIVER     1
setenv G4VIS_BUILD_DAWNFILE_DRIVER 1
setenv G4VIS_BUILD_VRML_DRIVER     1
setenv G4VIS_BUILD_VRMLFILE_DRIVER 1
# G4 use flags :
setenv G4UI_USE_XM                 1
setenv G4VIS_USE_OPENGLXM          1
setenv G4VIS_USE_OPENGLX           1
setenv G4VIS_USE_RAYX              1
setenv G4VIS_USE_OIX               1
# Specific :
setenv CLHEP_BASE_DIR /geant4/HP-UX
setenv OGLHOME        /geant4/HP-UX
setenv OIVFLAGS       "-I/geant4/HP-UX/include/SoFree"
setenv OIVLIBS        "-L/geant4/HP-UX/lib -lhepvisXt -lhepvis -lSoFreeXt -lSoFree"
setenv SOFREEUSER     /projects/SoFree/user/
# OPACS :
#setenv G4UI_BUILD_WO_SESSION       1
#setenv G4VIS_BUILD_OPACS_DRIVER    1
#setenv G4UI_USE_WO                 1
#setenv G4VIS_USE_OPACS             1
#setenv OCONFIG HP-UX-aCC
#source /projects/OPACS/setup.csh
endif
#---------------------------------------------------
if ( `uname -n` == asc ) then
setenv CVSROOT :pserver:barrand@g4cvs.cern.ch:/afs/cern.ch/sw/geant4/cvs
setenv G4INSTALL                   /geant4/stt/$REF/src/geant4
setenv G4SYSTEM                    OSF1
setenv G4WORKDIR                   /geant4/stt/$REF/$G4SYSTEM/$DEBOPT
setenv G4DEBUG                     1
#setenv G4MAKESHLIB                 $G4INSTALL/config/makeshlib.sh
# G4 build flags :
setenv G4UI_BUILD_XM_SESSION       1
setenv G4UI_BUILD_XAW_SESSION      1
setenv G4VIS_BUILD_OPENGLXM_DRIVER 1
setenv G4VIS_BUILD_OPENGLX_DRIVER  1
setenv G4VIS_BUILD_RAYX_DRIVER     1
setenv G4VIS_BUILD_OIX_DRIVER      1
setenv G4VIS_BUILD_DAWN_DRIVER     1
setenv G4VIS_BUILD_DAWNFILE_DRIVER 1
setenv G4VIS_BUILD_VRML_DRIVER     1
setenv G4VIS_BUILD_VRMLFILE_DRIVER 1
# G4 use flags :
setenv G4UI_USE_XM                 1
setenv G4UI_USE_XAW                1
setenv G4VIS_USE_OPENGLXM          1
setenv G4VIS_USE_OPENGLX           1
setenv G4VIS_USE_RAYX              1
setenv G4VIS_USE_OIX               1
# Specific :
setenv RWBASE         /geant4/OSF1
setenv CLHEP_BASE_DIR /geant4/OSF1
setenv OGLHOME        /geant4/OSF1
setenv OIVHOME        /geant4/OpenInventor2.4.1
setenv OIVFLAGS       "-I$OIVHOME/include -I/geant4/OSF1/include"
setenv OIVLIBS        "-L/geant4/OSF1/lib -lhepvisXt -lhepvis  -L$OIVHOME/lib -lInventorXt -lInventor -limage"
setenv XENVIRONMENT   $OIVHOME/app-defaults/Inventor
# OPACS :
#setenv G4UI_BUILD_WO_SESSION       1
#setenv G4VIS_BUILD_OPACS_DRIVER    1
#setenv G4UI_USE_WO                 1
#setenv G4VIS_USE_OPACS             1
#source /projects/OPACS/setup.csh
endif
####################################################
####################################################
