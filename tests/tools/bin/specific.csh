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
# Steve O'Neale Birmingham/Cern ( CLHEP dev )
# Steve O'Neale Birmingham/Cern ( Reactivate persistancy)
####################################################

setenv REF undefined
if ( `pwd | grep /stt/dev1/` != "" ) then
  setenv REF dev1
endif
if ( `pwd | grep /stt/dev2/` != "" ) then
  setenv REF dev2
endif
if ( `pwd | grep /stt/prod/` != "" ) then
  setenv REF prod
endif

if ( $?G4DEBUG ) then
  setenv DEBOPT debug
else
  setenv DEBOPT optim
endif

# Generaal G4 build flags :
setenv G4UI_BUILD_TERMINAL_SESSION 1
setenv G4UI_BUILD_GAG_SESSION      1
setenv G4VIS_BUILD_DAWN_DRIVER     1
setenv G4VIS_BUILD_DAWNFILE_DRIVER 1
setenv G4VIS_BUILD_VRML_DRIVER     1
setenv G4VIS_BUILD_VRMLFILE_DRIVER 1
setenv G4VIS_BUILD_RAYTRACER_DRIVER 1

if ( `uname -n | grep rsplus` != "" ) then
  setenv G4USE_OSPACE 1
  setenv CVSROOT /afs/cern.ch/sw/geant4/cvs
  setenv G4SYSTEM AIX-xlC
  setenv G4INSTALL /afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  setenv G4STTDIR  /afs/cern.ch/sw/geant4/stt/$REF/testtools/geant4/tests/tools
  setenv G4WORKDIR  /afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  setenv G4LIB $G4WORKDIR/lib
  # Other G4 build flags :
  #####setenv G4UI_BUILD_XM_SESSION       1
  #####setenv G4VIS_BUILD_OPENGLXM_DRIVER 1
  setenv G4VIS_BUILD_OPENGLX_DRIVER  1
  setenv XKEYSYMDB /usr/lib/X11/XKeysymDB
  setenv OGLHOME /afs/cern.ch/sw/geant4/dev/Mesa/Mesa-1.2.8
  ##### setenv G4VIS_BUILD_OIX_DRIVER      1
endif

if ( `uname -n | grep dxplus` != "" || \
     `uname -n | grep dcosf01` != "" ) then
  if ( $?G4STTNONISO ) then
    setenv DEBOPT ${DEBOPT}_NONISO
    setenv G4USE_OSPACE 1
    setenv CLHEP_BASE_DIR /afs/cern.ch/sw/geant4/dev/CLHEP/DEC-cxx/pro
  else
    setenv DEBOPT ${DEBOPT}_ISO
    setenv CLHEP_BASE_DIR /afs/cern.ch/sw/geant4/dev/CLHEP/DEC-cxx/iso
  endif
  if ( `uname -n | grep dcosf01` != "" ) then
    setenv CLHEP_LIB CLHEP-cxx62
  endif
  setenv CVSROOT /afs/cern.ch/sw/geant4/cvs
  setenv G4SYSTEM DEC-cxx
  setenv G4INSTALL /afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  setenv G4STTDIR  /afs/cern.ch/sw/geant4/stt/$REF/testtools/geant4/tests/tools
  setenv G4WORKDIR  /afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  setenv G4LIB $G4WORKDIR/lib
  # G4 build flags :
  #####setenv G4UI_BUILD_XM_SESSION       1
  #####setenv G4VIS_BUILD_OPENGLXM_DRIVER 1
  setenv G4VIS_BUILD_OPENGLX_DRIVER  1
  setenv XKEYSYMDB /usr/lib/X11/XKeysymDB
  setenv OGLHOME /afs/cern.ch/sw/geant4/dev/Mesa/Mesa-1.2.8
  ##### setenv G4VIS_BUILD_OIX_DRIVER      1
endif

if ( `uname -n | grep pcgeant`   != "" || \
     `uname -n | grep pcg4speed` != "" ) then
  setenv CVSROOT /afs/cern.ch/sw/geant4/cvs
  setenv G4SYSTEM Linux-g++
  setenv G4INSTALL /afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  setenv G4STTDIR  /afs/cern.ch/sw/geant4/stt/$REF/testtools/geant4/tests/tools
  setenv G4WORKDIR  /afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  setenv G4LIB $G4WORKDIR/lib
  # G4 build flags :
  #####setenv G4UI_BUILD_XM_SESSION       1
  #####setenv G4VIS_BUILD_OPENGLXM_DRIVER 1
  setenv G4VIS_BUILD_OPENGLX_DRIVER  1
  setenv XKEYSYMDB /usr/lib/X11/XKeysymDB
  setenv OGLHOME /afs/cern.ch/sw/geant4/dev/Mesa/Mesa-1.2.8
  ##### setenv G4VIS_BUILD_OIX_DRIVER      1
   setenv CLHEP_BASE_DIR /afs/cern.ch/sw/geant4/dev/CLHEP/$G4SYSTEM/pro
endif

if ( `uname -n | grep sgmedia` != "" ) then
  setenv G4USE_OSPACE 1
  setenv CVSROOT /afs/cern.ch/sw/geant4/cvs
  setenv G4SYSTEM SGI-CC
  setenv G4INSTALL /afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  setenv G4STTDIR  /afs/cern.ch/sw/geant4/stt/$REF/testtools/geant4/tests/tools
  setenv G4WORKDIR  /afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  setenv G4LIB $G4WORKDIR/lib
  # G4 build flags :
  #######setenv G4UI_BUILD_XM_SESSION       1
  setenv G4VIS_BUILD_OPENGLXM_DRIVER 1
  setenv G4VIS_BUILD_OPENGLX_DRIVER  1
  setenv G4VIS_BUILD_OIX_DRIVER      1
endif

if ( `uname -n | grep sun` != "" ) then
  setenv CVSROOT /afs/cern.ch/sw/geant4/cvs
  setenv G4INSTALL /afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  setenv G4STTDIR  /afs/cern.ch/sw/geant4/stt/$REF/testtools/geant4/tests/tools
  if ( $?G4STTNONISO ) then
    setenv G4SYSTEM SUN-CC
    setenv DEBOPT ${DEBOPT}_NONISO
    setenv G4USE_OSPACE 1
    setenv PATH `echo $PATH | sed s/SUNWspro50/SUNWspro/`
    setenv CLHEP_BASE_DIR /afs/cern.ch/sw/geant4/dev/CLHEP/$G4SYSTEM/pro
    # Persistency...
    if ( $?G4USE_HEPODBMS ) then  # Protect against double calling.
    else
      source $G4INSTALL/examples/extended/persistency/PersistentEx01/g4odbms_setup.csh
      setenv G4EXAMPLE_FDID 207
    endif
  else
    setenv G4SYSTEM SUN-CC5
    setenv DEBOPT ${DEBOPT}_ISO
    unsetenv G4USE_OSPACE
    setenv PATH `echo $PATH | sed s/SUNWspro50/SUNWspro/`
    setenv PATH `echo $PATH | sed s/SUNWspro/SUNWspro50/`
    # Persistency...
    unsetenv G4USE_HEPODBMS
    setenv CLHEP_BASE_DIR /afs/cern.ch/sw/geant4/dev/CLHEP/$G4SYSTEM/pro
  endif
  setenv G4WORKDIR  /afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  setenv G4LIB $G4WORKDIR/lib
  # G4 build flags :
  #######setenv G4UI_BUILD_XM_SESSION       1
  #######setenv G4VIS_BUILD_OPENGLXM_DRIVER 1
  #######setenv G4VIS_BUILD_OPENGLX_DRIVER  1
  #######setenv G4VIS_BUILD_OIX_DRIVER      1
endif

if ( `uname -n | grep hp` != "" ) then
  setenv G4USE_OSPACE 1
  setenv CVSROOT /afs/cern.ch/sw/geant4/cvs
  setenv G4SYSTEM HP-aCC
  setenv G4INSTALL /afs/cern.ch/sw/geant4/stt/$REF/src/geant4
  setenv G4STTDIR  /afs/cern.ch/sw/geant4/stt/$REF/testtools/geant4/tests/tools
  setenv G4WORKDIR  /afs/cern.ch/sw/geant4/stt/$REF/$G4SYSTEM/$DEBOPT
  setenv G4LIB $G4WORKDIR/lib
  setenv CLHEP_BASE_DIR /afs/cern.ch/sw/geant4/dev/CLHEP/$G4SYSTEM/pro
  # G4 build flags :
  ######setenv G4UI_BUILD_XM_SESSION       1
  setenv G4VIS_BUILD_OPENGLXM_DRIVER 1
  setenv G4VIS_BUILD_OPENGLX_DRIVER  1
  setenv G4VIS_BUILD_OIX_DRIVER      1
endif

####################################################
### LAL Orsay machines                           ###
####################################################
if ( `uname -n` == aleph ) then
  setenv CVSROOT :pserver:barrand@g4cvs.cern.ch:/afs/cern.ch/sw/geant4/cvs
  setenv G4INSTALL /geant4/geant4-01-00-ref-02
  setenv G4WORKDIR $G4INSTALL
  setenv G4STTDIR  $G4WORKDIR/stt
  setenv G4LIB     $G4WORKDIR/lib
  setenv G4SYSTEM  HP-aCC
  setenv G4DEBUG   1
  #setenv G4MAKESHLIB                 $G4INSTALL/config/makeshlib.sh
  # G4 build flags :
  setenv G4UI_BUILD_XM_SESSION       1
  setenv G4VIS_BUILD_OPENGLXM_DRIVER 1
  setenv G4VIS_BUILD_OPENGLX_DRIVER  1
  setenv G4VIS_BUILD_OIX_DRIVER      1
  setenv G4VIS_BUILD_DAWN_DRIVER     1
  setenv G4VIS_BUILD_DAWNFILE_DRIVER 1
  setenv G4VIS_BUILD_VRML_DRIVER     1
  setenv G4VIS_BUILD_VRMLFILE_DRIVER 1
  # G4 use flags :
  setenv G4UI_USE_XM                 1
  setenv G4VIS_USE_OPENGLXM          1
  setenv G4VIS_USE_OPENGLX           1
  setenv G4VIS_USE_OIX               1
  setenv G4VIS_USE_DAWN              1
  setenv G4VIS_USE_DAWNFILE          1
  setenv G4VIS_USE_VRML              1
  setenv G4VIS_USE_VRMLFILE          1
  # Specific :
  setenv CLHEP_BASE_DIR /lal/CLHEP/1.4/HP-UX-aCC
  setenv OGLHOME        /lal/Mesa/3.1/HP-UX
  setenv OGLLIBS        "-L$OGLHOME/lib -lGLU -lGL"
  setenv OIVHOME        /lal/SoFree/v2r9
  setenv HEPVISHOME     /lal/HEPVis/v5r1-05-LAL
  setenv OIVFLAGS       "-I$OIVHOME/include -I$HEPVISHOME/include"
  setenv OIVLIBS        "-L$HEPVISHOME/HP-UX-aCC-SF -lHEPVis -L$OIVHOME/HP-UX-aCC -lSoFree"
  setenv SOFREEUSER     $OIVHOME/user/
  # OPACS :
  setenv G4UI_BUILD_WO_SESSION       1
  setenv G4VIS_BUILD_OPACS_DRIVER    1
  setenv G4UI_USE_WO                 1
  setenv G4VIS_USE_OPACS             1
  setenv OCONFIG HP-UX
  source /lal/OPACS/v3/setup.csh
  # Else :
  #setenv XENVIRONMENT   $G4INSTALL/tests/test201/test201.xrm
  setenv XENVIRONMENT   g4.xrm
  setenv PATH "${PATH}:/lal/DAWN/3.72b/HP-UX-aCC"
  setenv CPPVERBOSE 1
  set prompt='g4-aleph> ' 
endif
#---------------------------------------------------
if ( `uname -n` == asc ) then
  # In CLHEP-default.h : //GB #define HEP_USE_STD 1
  #                      //GB #define HEP_HAVE_BOOL 1
  setenv CVSROOT :pserver:barrand@g4cvs.cern.ch:/afs/cern.ch/rd44/cvs
  setenv G4INSTALL /geant4/geant4-02-00
  setenv G4WORKDIR $G4INSTALL
  setenv G4STTDIR  $G4WORKDIR/stt
  setenv G4LIB     $G4WORKDIR/lib
  setenv G4SYSTEM DEC-cxx
  setenv G4DEBUG 1
  #setenv G4MAKESHLIB $G4INSTALL/config/makeshlib.sh
  # G4 build flags :
  setenv G4UI_BUILD_XM_SESSION       1
  setenv G4UI_BUILD_XAW_SESSION      1
  setenv G4VIS_BUILD_OPENGLXM_DRIVER 1
  setenv G4VIS_BUILD_OPENGLX_DRIVER  1
  #TGS out with cxx ansi : 
  #setenv G4VIS_BUILD_OIX_DRIVER      1
  setenv G4VIS_BUILD_DAWN_DRIVER     1
  setenv G4VIS_BUILD_DAWNFILE_DRIVER 1
  # Pb with cxx ansi :
  #setenv G4VIS_BUILD_VRML_DRIVER     1
  #setenv G4VIS_BUILD_VRMLFILE_DRIVER 1
  # G4 use flags :
  setenv G4UI_USE_XM                 1
  setenv G4UI_USE_XAW                1
  setenv G4VIS_USE_OPENGLXM          1
  setenv G4VIS_USE_OPENGLX           1
  #setenv G4VIS_USE_OIX               1
  setenv G4VIS_USE_DAWN              1
  setenv G4VIS_USE_DAWNFILE          1
  # Pb with cxx ansi :
  #setenv G4VIS_USE_VRML              1
  #setenv G4VIS_USE_VRMLFILE          1
  # Specific :
  setenv CLHEP_BASE_DIR /lal/CLHEP/1.5/OSF1-cxx
  setenv OGLHOME        /lal/Mesa/3.2/OSF1
  setenv OIVHOME        /lal/OpenInventor/2.5
  setenv HEPVISHOME     /lal/HEPVis/v5r1-05-LAL
  setenv OIVFLAGS       "-I$OIVHOME/include -I$HEPVISHOME/include"
  setenv OIVLIBS        "-L$HEPVISHOME/OSF1-cxx-TGS -lHEPVis -L$OIVHOME/lib -lInventorXt -lInventor -limage"
  #setenv XENVIRONMENT   $OIVHOME/app-defaults/Inventor
  setenv XENVIRONMENT   g4.xrm
  # OPACS :
  setenv G4UI_BUILD_WO_SESSION       1
  setenv G4VIS_BUILD_OPACS_DRIVER    1
  setenv G4UI_USE_WO                 1
  setenv G4VIS_USE_OPACS             1
  source /lal/OPACS/v3/setup.csh
  set prompt=${G4INSTALL}-${G4SYSTEM}'> '
  # Else :
  setenv LD_LIBRARY_PATH "$OIVHOME/lib:$HEPVISHOME/OSF1-cxx-TGS:${OGLHOME}/lib"
  setenv PATH "${PATH}:/lal/DAWN/3.72b/OSF1-cxx"
  setenv XENVIRONMENT   g4Xt.xrm
  setenv CPPVERBOSE 1
  set prompt='g4-asc> ' 
  # To be able to link :
  limit datasize 500000
endif
#---------------------------------------------------
if ( `uname -n` == "lx1.lal.in2p3.fr" ) then
  setenv CVSROOT :pserver:barrand@g4cvs.cern.ch:/afs/cern.ch/sw/geant4/cvs
  setenv G4INSTALL /geant4/geant4-01-01
  setenv G4WORKDIR $G4INSTALL
  setenv G4STTDIR  $G4WORKDIR/stt
  setenv G4LIB     $G4WORKDIR/lib
  setenv G4SYSTEM  Linux-g++
  setenv G4DEBUG   1
  #setenv G4MAKESHLIB                 $G4INSTALL/config/makeshlib.sh
  # G4 build flags :
  setenv G4UI_BUILD_XM_SESSION       1
  setenv G4VIS_BUILD_OPENGLXM_DRIVER 1
  setenv G4VIS_BUILD_OPENGLX_DRIVER  1
  setenv G4VIS_BUILD_OIX_DRIVER      1
  setenv G4VIS_BUILD_DAWN_DRIVER     1
  setenv G4VIS_BUILD_DAWNFILE_DRIVER 1
  setenv G4VIS_BUILD_VRML_DRIVER     1
  setenv G4VIS_BUILD_VRMLFILE_DRIVER 1
  # G4 use flags :
  setenv G4UI_USE_XM                 1
  setenv G4VIS_USE_OPENGLXM          1
  setenv G4VIS_USE_OPENGLX           1
  setenv G4VIS_USE_OIX               1
  setenv G4VIS_USE_DAWN              1
  setenv G4VIS_USE_DAWNFILE          1
  setenv G4VIS_USE_VRML              1
  setenv G4VIS_USE_VRMLFILE          1
  # Specific :
  setenv CLHEP_BASE_DIR /lal/CLHEP/1.4/Linux-gxx
  setenv OGLHOME        /lal/Mesa/3.1/Linux
  setenv OIVHOME        /lal/SoFree/v2r9
  setenv HEPVISHOME     /lal/HEPVis/v5r1-05-LAL
  setenv OIVFLAGS       "-I$OIVHOME/include -I$HEPVISHOME/include"
  setenv OIVLIBS        "-L$HEPVISHOME/Linux-gxx-SF -lHEPVis -L$OIVHOME/Linux-gxx -lSoFree"
  setenv SOFREEUSER     $OIVHOME/user/
  # OPACS :
  setenv G4UI_BUILD_WO_SESSION       1
  setenv G4VIS_BUILD_OPACS_DRIVER    1
  setenv G4UI_USE_WO                 1
  setenv G4VIS_USE_OPACS             1
  source /lal/OPACS/v3/setup.csh
  setenv LD_LIBRARY_PATH "${LD_LIBRARY_PATH}:$OIVHOME/Linux-gxx:$HEPVISHOME/Linux-gxx-SF:${OGLHOME}/lib"
  # Else :
  #setenv XENVIRONMENT   $G4INSTALL/tests/test201/test201.xrm
  setenv XENVIRONMENT   g4.xrm
  setenv PATH "${PATH}:/lal/DAWN/3.72b/Linux-egcs"
  setenv CPPVERBOSE 1
  set prompt='g4-lx1> ' 
endif
#
if ( `uname -n` == "papou1" ) then
  setenv CVSROOT :pserver:barrand@g4cvs.cern.ch:/afs/cern.ch/sw/geant4/cvs
  setenv G4INSTALL /geant4/geant4-01-00-ref-02
  setenv G4WORKDIR $G4INSTALL
  setenv G4STTDIR  $G4WORKDIR/stt
  setenv G4LIB     $G4WORKDIR/lib
  setenv G4SYSTEM  SUN-CC5
  setenv G4DEBUG   1
  #setenv G4MAKESHLIB                 $G4INSTALL/config/makeshlib.sh
  # G4 build flags :
  setenv G4UI_BUILD_XM_SESSION       1
  setenv G4VIS_BUILD_OPENGLXM_DRIVER 1
  setenv G4VIS_BUILD_OPENGLX_DRIVER  1
  #setenv G4VIS_BUILD_OIX_DRIVER      1
  setenv G4VIS_BUILD_DAWN_DRIVER     1
  setenv G4VIS_BUILD_DAWNFILE_DRIVER 1
  setenv G4VIS_BUILD_VRML_DRIVER     1
  setenv G4VIS_BUILD_VRMLFILE_DRIVER 1
  # G4 use flags :
  setenv G4UI_USE_XM                 1
  setenv G4VIS_USE_OPENGLXM          1
  setenv G4VIS_USE_OPENGLX           1
  #setenv G4VIS_USE_OIX               1
  # Specific :
  setenv CLHEP_BASE_DIR /lal/CLHEP/1.4/SunOS-CC
  setenv OGLHOME        /lal/Mesa/3.0/SunOS
  setenv OIVHOME        /lal/SoFree/v2r9
  setenv HEPVISHOME     /lal/HEPVis/v5r1-05-LAL
  setenv OIVFLAGS       "-I$OIVHOME/include -I$HEPVISHOME/include"
  setenv OIVLIBS        "-L$HEPVISHOME/SunOS-CC -lHEPVis -L$OIVHOME/SunOS-CC -lSoFree"
  setenv SOFREEUSER     $OIVHOME/user/
  # OPACS :
  #setenv G4UI_BUILD_WO_SESSION       1
  #setenv G4VIS_BUILD_OPACS_DRIVER    1
  #setenv G4UI_USE_WO                 1
  #setenv G4VIS_USE_OPACS             1
  #source /lal/OPACS/v3/setup.csh
  #setenv LD_LIBRARY_PATH "$OIVHOME/Linux-gxx:/lal/HEPVis/v5r0/Linux-gxx"
  # Else :
  setenv CPPVERBOSE 1
  set prompt='g4-papou1> ' 
endif
# Comments :
# --------
#
if ( `uname -n` == "pc-gbp" || `uname -n` == "pc100" ) then
  setenv CVSROOT :pserver:barrand@g4cvs.cern.ch:/afs/cern.ch/sw/geant4/cvs
  setenv G4INSTALL /geant4/geant4-02-00
  setenv G4WORKDIR $G4INSTALL
  setenv G4STTDIR  $G4WORKDIR/stt
  setenv G4LIB     $G4WORKDIR/lib
  setenv G4SYSTEM  Linux-g++
  setenv G4DEBUG   1
  #setenv G4MAKESHLIB                 $G4INSTALL/config/makeshlib.sh
  # G4 build flags :
  setenv G4UI_BUILD_XM_SESSION       1
  setenv G4VIS_BUILD_OPENGLXM_DRIVER 1
  setenv G4VIS_BUILD_OPENGLX_DRIVER  1
  setenv G4VIS_BUILD_OIX_DRIVER      1
  setenv G4VIS_BUILD_DAWN_DRIVER     1
  setenv G4VIS_BUILD_DAWNFILE_DRIVER 1
  setenv G4VIS_BUILD_VRML_DRIVER     1
  setenv G4VIS_BUILD_VRMLFILE_DRIVER 1
  setenv G4ANALYSIS_BUILD_OPEN_SCIENTIST 1
  setenv G4ANALYSIS_BUILD_JAS            1
  # G4 use flags :
  setenv G4UI_USE_XM                 1
  setenv G4VIS_USE_OPENGLXM          1
  setenv G4VIS_USE_OPENGLX           1
  setenv G4VIS_USE_OIX               1
  setenv G4VIS_USE_DAWN              1
  setenv G4VIS_USE_DAWNFILE          1
  setenv G4VIS_USE_VRML              1
  setenv G4VIS_USE_VRMLFILE          1
  setenv G4ANALYSIS_USE_OPEN_SCIENTIST 1
  setenv G4ANALYSIS_USE_JAS            1
  # Specific :
  setenv CLHEP_BASE_DIR /lal/CLHEP/1.5/Linux-gxx
  setenv OGLHOME        /usr/X11R6
  #  setenv OGLHOME        /lal/Mesa/3.2/Linux
  setenv OIVHOME        /lal/SoFree/v2r9
  setenv HEPVISHOME     /lal/HEPVis/v5r1-05-LAL
  setenv OIVFLAGS       "-I$OIVHOME/include -I$HEPVISHOME/include"
  setenv OIVLIBS        "-L$HEPVISHOME/Linux-gxx-SF -lHEPVis -L$OIVHOME/Linux-gxx -lSoFree"
  setenv SOFREEUSER     $OIVHOME/user/
  # OPACS :
  setenv G4UI_BUILD_WO_SESSION       1
  setenv G4VIS_BUILD_OPACS_DRIVER    1
  setenv G4UI_USE_WO                 1
  setenv G4VIS_USE_OPACS             1
  source /lal/OPACS/v3/setup.csh
  setenv JDKHOME /lal/JDK/jdk1.2.2
  setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:$OIVHOME/Linux-gxx:$HEPVISHOME/Linux-gxx-SF:${OGLHOME}/lib
  setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:$JDKHOME/jre/lib/i386:$JDKHOME/jre/lib/i386/classic:$JDKHOME/jre/lib/i386/native_threads
set jars=/lal/jas/2.0alpha4/release/lib
  setenv CLASSPATH ${CLASSPATH}:$jars/collections.jar:$jars/hep.jar:$jars/jas.jar
  # OpenScientist :
  source /projects/OpenScientist/v6r0/omake/setup.csh
  setenv AIDAINCS $HCLHOME/include
  # Else :
  setenv XENVIRONMENT   g4Xt.xrm
  setenv PATH "${PATH}:/lal/DAWN/dawn_3_85a/Linux/bin"
  setenv CPPVERBOSE 1
  alias g4ANA01 "cd $G4INSTALL/examples/extended/analysis/AnaEx01"
  alias ana01   "$G4WORKDIR/bin/$G4SYSTEM/AnaEx01"
  set prompt='g4-pc-gbp> ' 
endif
#
####################################################
####################################################
