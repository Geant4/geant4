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

if ( $?G4DEBUG ) then
  setenv DEBOPT debug
else
  setenv DEBOPT optim
endif

if ( $?G4USE_STL ) then
  setenv DEBOPT ${DEBOPT}_STL
  if ( `uname -n | grep sun` != "" ) then
    setenv G4USE_OSPACE 1
  endif
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

if ( `uname -n | grep pcgeant` != "" ) then
  setenv CVSROOT /afs/cern.ch/sw/geant4/cvs
  setenv G4SYSTEM Linux-g++
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
  setenv LD_LIBRARY_PATH /afs/cern.ch/sw/geant4/dev/CLHEP/Linux-g++/new/lib
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
  setenv G4VIS_BUILD_OPENGLXM_DRIVER 1
  setenv G4VIS_BUILD_OPENGLX_DRIVER  1
  setenv G4VIS_BUILD_RAYX_DRIVER     1
  setenv G4VIS_BUILD_OIX_DRIVER      1
  setenv G4VIS_BUILD_DAWN_DRIVER     1
  setenv G4VIS_BUILD_DAWNFILE_DRIVER 1
  setenv G4VIS_BUILD_VRML_DRIVER     1
  setenv G4VIS_BUILD_VRMLFILE_DRIVER 1
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
  setenv G4VIS_BUILD_DAWN_DRIVER     1
  setenv G4VIS_BUILD_DAWNFILE_DRIVER 1
  setenv G4VIS_BUILD_VRML_DRIVER     1
  setenv G4VIS_BUILD_VRMLFILE_DRIVER 1
  setenv LD_LIBRARY_PATH /opt/SUNWspro/lib:/usr/openwin/lib:/usr/dt/lib:/usr/local/lib:/afs/cern.ch/sw/geant4/dev/CLHEP/SUN-CC/new/lib
endif

if ( `uname -n | grep hp` != "" ) then
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
setenv CVSROOT :pserver:barrand@g4cvs.cern.ch:/afs/cern.ch/rd44/cvs
setenv G4INSTALL /geant4/dev/geant4
setenv G4LIB     /geant4/dev/lib
setenv G4WORKDIR /geant4/dev
setenv G4SYSTEM  HP-aCC
setenv G4DEBUG   1
#setenv G4MAKESHLIB                 $G4INSTALL/config/makeshlib.sh
# G4 build flags :
setenv G4UI_BUILD_XM_SESSION       1
setenv G4VIS_BUILD_OPENGLXM_DRIVER 1
setenv G4VIS_BUILD_OPENGLX_DRIVER  1
setenv G4VIS_BUILD_RAYX_DRIVER     1
setenv G4VIS_BUILD_OIX_DRIVER      1
#setenv G4VIS_BUILD_DAWN_DRIVER     1
#setenv G4VIS_BUILD_DAWNFILE_DRIVER 1
#setenv G4VIS_BUILD_VRML_DRIVER     1
#setenv G4VIS_BUILD_VRMLFILE_DRIVER 1
# G4 use flags :
setenv G4UI_USE_XM                 1
setenv G4VIS_USE_OPENGLXM          1
setenv G4VIS_USE_OPENGLX           1
setenv G4VIS_USE_RAYX              1
setenv G4VIS_USE_OIX               1
# Specific :
setenv CLHEP_BASE_DIR /lal/CLHEP/1.4/HP-UX-aCC
setenv OGLHOME        /lal/Mesa/3.0/HP-UX
setenv OGLLIBS        "-L$OGLHOME/lib -lGLU -lGL"
#setenv OIVHOME        /lal/SoFree/v2r3
#setenv HEPVISHOME     /lal/HEPVis/v5r0-06-LAL
setenv OIVHOME        /projects/SoFree
setenv HEPVISHOME     /projects/HEPVis
setenv OIVFLAGS       "-I$OIVHOME/include -I$HEPVISHOME/include"
setenv OIVLIBS        "-L$HEPVISHOME/HP-UX-aCC-SF -lHEPVis -L$OIVHOME/HP-UX-aCC -lSoFree"
setenv SOFREEUSER     $OIVHOME/user/
# OPACS :
setenv G4UI_BUILD_WO_SESSION       1
setenv G4VIS_BUILD_OPACS_DRIVER    1
setenv G4UI_USE_WO                 1
setenv G4VIS_USE_OPACS             1
setenv OCONFIG HP-UX-aCC
source /lal/OPACS/v3/setup.csh
# Else :
setenv CPPVERBOSE 1
endif
#---------------------------------------------------
if ( `uname -n` == asc ) then
# In CLHEP-default.h : //#define HEP_USE_STD 1
setenv CVSROOT :pserver:barrand@g4cvs.cern.ch:/afs/cern.ch/rd44/cvs
setenv G4INSTALL /geant4/dev/geant4
setenv G4LIB     /geant4/dev/lib
setenv G4WORKDIR /geant4/dev
setenv G4SYSTEM DEC-cxx
setenv G4DEBUG 1
#setenv G4MAKESHLIB $G4INSTALL/config/makeshlib.sh
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
setenv RWBASE         /lal/rogue/6.1/OSF1-cxx
setenv CLHEP_BASE_DIR /lal/CLHEP/1.4/OSF1-cxx
setenv OGLHOME        /lal/Mesa/3.0/OSF1
setenv OIVHOME        /lal/OpenInventor/2.5
setenv HEPVISHOME     /projects/HEPVis
setenv OIVFLAGS       "-I$OIVHOME/include -I$HEPVISHOME/include"
setenv OIVLIBS        "-L$HEPVISHOME/OSF1-cxx-TGS -lHEPVis -L$OIVHOME/lib -lInventorXt -lInventor -limage"
setenv XENVIRONMENT   $OIVHOME/app-defaults/Inventor
# OPACS :
setenv G4UI_BUILD_WO_SESSION       1
setenv G4VIS_BUILD_OPACS_DRIVER    1
setenv G4UI_USE_WO                 1
setenv G4VIS_USE_OPACS             1
source /lal/OPACS/v3/setup.csh
set prompt=${G4INSTALL}-${G4SYSTEM}'> '
# Else :
setenv CPPVERBOSE 1
endif
#---------------------------------------------------
if ( `uname -n` == "lx1.lal.in2p3.fr" ) then
setenv CVSROOT :pserver:barrand@g4cvs.cern.ch:/afs/cern.ch/rd44/cvs
setenv G4INSTALL /geant4/dev/geant4
setenv G4LIB     /geant4/dev/lib
setenv G4WORKDIR /geant4/dev
setenv G4SYSTEM  Linux-g++
setenv G4DEBUG   1
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
setenv RWBASE         /lal/rogue/6.1/Linux-gxx
setenv CLHEP_BASE_DIR /lal/CLHEP/1.4/Linux-gxx
setenv OGLHOME        /lal/Mesa/3.0/Linux
#setenv OIVHOME        /lal/SoFree/v2r3
#setenv HEPVISHOME     /lal/HEPVis/v5r0-06-LAL
setenv OIVHOME        /projects/SoFree
setenv HEPVISHOME     /projects/HEPVis
setenv OIVFLAGS       "-I$OIVHOME/include -I$HEPVISHOME/include"
setenv OIVLIBS        "-L$HEPVISHOME/Linux-egcs-SF -lHEPVis -L$OIVHOME/Linux-egcs -lSoFree"
setenv SOFREEUSER     $OIVHOME/user/
# OPACS :
setenv G4UI_BUILD_WO_SESSION       1
setenv G4VIS_BUILD_OPACS_DRIVER    1
setenv G4UI_USE_WO                 1
setenv G4VIS_USE_OPACS             1
setenv OCONFIG HP-UX-aCC
source /lal/OPACS/v3/setup.csh
setenv LD_LIBRARY_PATH "${LD_LIBRARY_PATH}:$OIVHOME/Linux-egcs:/lal/HEPVis/v5r0-06-LAL/Linux-egcs-SF"
# Else :
setenv CPPVERBOSE 1
endif
# Comments :
# --------
#  Problems with advance method in rogue/6.1/rw/tpsldict.h, 
# tvsldict.h. Solved by editing directly the file !
#---------------------------------------------------
if ( `uname -n` == "papou1" ) then
setenv CVSROOT :pserver:barrand@g4cvs.cern.ch:/afs/cern.ch/rd44/cvs
setenv G4INSTALL /geant4/dev/geant4
setenv G4LIB     /geant4/dev/lib
setenv G4WORKDIR /geant4/dev
setenv G4SYSTEM  SUN-CC
setenv G4DEBUG   1
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
setenv RWBASE         /lal/rogue/6.1/SunOS-CC
setenv CLHEP_BASE_DIR /lal/CLHEP/1.4/SunOS-CC
setenv OGLHOME        /lal/Mesa/3.0/SunOS
setenv OIVHOME        /lal/SoFree/v2r3
setenv OIVFLAGS       "-I$OIVHOME/include -I/lal/HEPVis/v5r0-06-LAL/include"
setenv OIVLIBS        "-L/lal/HEPVis/v5r0-06-LAL/SunOS-CC-SF -lHEPVis -L$OIVHOME/SunOS-CC -lSoFree"
setenv SOFREEUSER     $OIVHOME/user/
# OPACS :
setenv G4UI_BUILD_WO_SESSION       1
setenv G4VIS_BUILD_OPACS_DRIVER    1
setenv G4UI_USE_WO                 1
setenv G4VIS_USE_OPACS             1
setenv OCONFIG HP-UX-aCC
source /lal/OPACS/v3/setup.csh
setenv LD_LIBRARY_PATH "$OIVHOME/Linux-gxx:/lal/HEPVis/v5r0/Linux-gxx"
# Else :
setenv CPPVERBOSE 1
endif
# Comments :
# --------
####################################################
####################################################
