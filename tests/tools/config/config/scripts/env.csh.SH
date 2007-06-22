case $CONFIG in
'')
	if test -f config.sh; then TOP=.;
	elif test -f ../config.sh; then TOP=..;
	elif test -f ../../config.sh; then TOP=../..;
	elif test -f ../../../config.sh; then TOP=../../..;
	elif test -f ../../../../config.sh; then TOP=../../../..;
	else
		echo "Can't find config.sh."; exit 1
	fi
	. $TOP/config.sh
	;;
esac
: This forces SH files to create target in same directory as SH file.
: This is so that make depend always knows where to find SH derivatives.
case "$0" in
*/*) cd `expr X$0 : 'X\(.*\)/'` ;;
esac
#echo "Extracting env.csh (with variable substitutions)"
: This section of the file will have variable substitutions done on it.
: Move anything that needs config subs from !NO!SUBS! section to !GROK!THIS!.
: Protect any dollar signs and backticks that you do not want interpreted
: by putting a backslash in front.  You may delete these comments.
$spitshell >env.csh <<!GROK!THIS!
######################################
#
# Clean all G4 envs
unsetenv  CLHEP_BASE_DIR
unsetenv  CLHEP_INCLUDE_DIR
unsetenv  CLHEP_LIB
unsetenv  CLHEP_LIB_DIR

unsetenv  G4ANALYSIS_USE
unsetenv  G4DEBUG
unsetenv  G4INCLUDE
unsetenv  G4INSTALL

unsetenv  G4LEDATA
unsetenv  G4LEVELGAMMADATA
unsetenv  G4NEUTRONHPDATA
unsetenv  G4RADIOACTIVEDATA

unsetenv  G4LIB
unsetenv  G4LIB_BUILD_G3TOG4
unsetenv  G4LIB_BUILD_SHARED
unsetenv  G4LIB_BUILD_STATIC
unsetenv  G4LIB_BUILD_ZLIB
unsetenv  G4LIB_USE_G3TOG4
unsetenv  G4LIB_USE_GRANULAR
unsetenv  G4LIB_USE_ZLIB

unsetenv  G4SYSTEM

unsetenv  G4UI_BUILD_WIN32_SESSION
unsetenv  G4UI_BUILD_XAW_SESSION
unsetenv  G4UI_BUILD_XM_SESSION
unsetenv  G4UI_USE_TCSH
unsetenv  G4UI_USE_WIN32
unsetenv  G4UI_USE_XAW
unsetenv  G4UI_USE_XM

unsetenv  G4VIS_BUILD_DAWN_DRIVER
unsetenv  G4VIS_BUILD_OIWIN32_DRIVER
unsetenv  G4VIS_BUILD_OIX_DRIVER
unsetenv  G4VIS_BUILD_OPENGLWIN32_DRIVER
unsetenv  G4VIS_BUILD_OPENGLXM_DRIVER
unsetenv  G4VIS_BUILD_OPENGLX_DRIVER
unsetenv  G4VIS_BUILD_RAYTRACERX_DRIVER
unsetenv  G4VIS_BUILD_VRML_DRIVER

unsetenv  G4VIS_USE_DAWN
unsetenv  G4VIS_USE_OIWIN32
unsetenv  G4VIS_USE_OIX
unsetenv  G4VIS_USE_OPENGLWIN32
unsetenv  G4VIS_USE_OPENGLX
unsetenv  G4VIS_USE_OPENGLXM
unsetenv  G4VIS_USE_RAYTRACERX
unsetenv  G4VIS_USE_VRML

######################################
#
# g4system.U
#
#+
setenv G4SYSTEM "$g4system"
echo "On this machine the G4SYSTEM=\$G4SYSTEM"

#
# g4dirs.U
#
#+
if ( X$g4install != X$g4final_install ) then
setenv G4INSTALL "$g4final_install/src/geant4"
echo "On this machine the G4INSTALL=\$G4INSTALL"
else
setenv G4INSTALL "$g4install"
echo "On this machine the G4INSTALL=\$G4INSTALL"
endif

#+
if ( X$g4includes_flag != Xn ) then
if ( X$g4install != X$g4final_install ) then
setenv G4INCLUDE "$g4final_install/include/geant4"
else
setenv G4INCLUDE "$g4final_install/include"
endif
echo "On this machine the G4INCLUDE=\$G4INCLUDE"
#else
#setenv G4BASE "$g4final_install/include/geant4/source"
#echo "On this machine the G4BASE=\$G4BASE"
endif

#+
if ( X$g4lib != X ) then
if ( X$g4install != X$g4final_install ) then
setenv G4LIB "$g4final_install/lib/geant4"
else
setenv G4LIB "$g4final_install/lib"
endif
echo "On this machine the G4LIB=\$G4LIB"
endif

#+
if ( X$g4levelgammadata != X ) then
setenv G4LEVELGAMMADATA "$g4levelgammadata"
echo "On this machine the G4LEVELGAMMADATA=\$G4LEVELGAMMADATA"
endif

#+
if ( X$g4radioactivedata != X ) then
setenv G4RADIOACTIVEDATA "$g4radioactivedata"
echo "On this machine the G4RADIOACTIVEDATA=\$G4RADIOACTIVEDATA"
endif

#+
if ( X$g4ledata != X ) then
setenv G4LEDATA "$g4ledata"
echo "On this machine the G4LEDATA=\$G4LEDATA"
endif

#+
if ( X$g4neutronhpcrosssections != X ) then
setenv G4NEUTRONHPDATA "$g4neutronhpcrosssections"
echo "On this machine the G4NEUTRONHPDATA=\$G4NEUTRONHPDATA"
endif

#+
#if ( X$g4elasticdata != X ) then
#setenv G4ELASTICDATA "$g4elasticdata"
#echo "On this machine the G4ELASTICDATA=\$G4ELASTICDATA"
#endif

#
# g4clhep.U
#
if ( X$g4clhep_base_dir != X ) then
setenv CLHEP_BASE_DIR "$g4clhep_base_dir"
echo "On this machine the CLHEP_BASE_DIR=\$CLHEP_BASE_DIR"
endif

#+
if ( X$g4clhep_include_dir != X ) then
setenv CLHEP_INCLUDE_DIR "$g4clhep_include_dir"
echo "On this machine the CLHEP_INCLUDE_DIR=\$CLHEP_INCLUDE_DIR"
endif

#+
if ( X$g4clhep_lib_dir != X ) then
setenv CLHEP_LIB_DIR "$g4clhep_lib_dir"
echo "On this machine the CLHEP_LIB_DIR=\$CLHEP_LIB_DIR"
endif

#+
if ( X$g4clhep_lib != X ) then
setenv CLHEP_LIB "$g4clhep_lib"
echo "On this machine the CLHEP_LIB=\$CLHEP_LIB"
endif

#+
#
# g4debug
#
if ( X$g4debug == Xy ) then
setenv G4DEBUG 1
echo "On this machine the G4DEBUG=\$G4DEBUG"
endif

#
# g4analysis
#
#+
if ( X$g4wanalysis_use == Xy ) then
setenv G4ANALYSIS_USE 1
echo "On this machine the G4ANALYSIS_USE=\$G4ANALYSIS_USE"
endif

#
# g4ui
#
#+
if ( X$g4ui_build_xaw_session == Xy ) then
setenv G4UI_BUILD_XAW_SESSION 1
echo "On this machine the G4UI_BUILD_XAW_SESSION=\$G4UI_BUILD_XAW_SESSION"
endif

#+
if ( X$g4ui_use_xaw == Xy ) then
setenv G4UI_USE_XAW 1
echo "On this machine the G4UI_USE_XAW=\$G4UI_USE_XAW"
endif

#+
if ( X$g4ui_build_xm_session == Xy ) then
setenv G4UI_BUILD_XM_SESSION 1
echo "On this machine the G4UI_BUILD_XM_SESSION=\$G4UI_BUILD_XM_SESSION"
endif

#+
if ( X$g4ui_use_xm == Xy ) then
setenv G4UI_USE_XM 1
echo "On this machine the G4UI_USE_XM=\$G4UI_USE_XM"
endif

#+
if ( X$g4ui_build_win32_session == Xy ) then
setenv G4UI_BUILD_WIN32_SESSION 1
echo "On this machine the G4UI_BUILD_WIN32_SESSION=\$G4UI_BUILD_WIN32_SESSION"
endif

#+
if ( X$g4ui_use_win32 == Xy ) then
setenv G4UI_USE_WIN32 1
echo "On this machine the G4UI_USE_WIN32=\$G4UI_USE_WIN32"
endif

#
# g4vis
#
#+
if ( X$g4vis_build_dawn_driver == Xy ) then
setenv G4VIS_BUILD_DAWN_DRIVER 1
echo "On this machine the G4VIS_BUILD_DAWN_DRIVER=\$G4VIS_BUILD_DAWN_DRIVER"
endif

#+
if ( X$g4vis_build_openglx_driver == Xy ) then
setenv G4VIS_BUILD_OPENGLX_DRIVER 1
echo "On this machine the G4VIS_BUILD_OPENGLX_DRIVER=\$G4VIS_BUILD_OPENGLX_DRIVER"
endif

#+
if ( X$g4vis_build_openglxm_driver == Xy ) then
setenv G4VIS_BUILD_OPENGLXM_DRIVER 1
echo "On this machine the G4VIS_BUILD_OPENGLXM_DRIVER=\$G4VIS_BUILD_OPENGLXM_DRIVER"
endif

#+
if ( X$g4vis_build_openglwin32_driver == Xy ) then
setenv G4VIS_BUILD_OPENGLWIN32_DRIVER 1
echo "On this machine the G4VIS_BUILD_OPENGLWIN32_DRIVER=\$G4VIS_BUILD_OPENGLWIN32_DRIVER"
endif

#+
if ( X$g4vis_build_oix_driver == Xy ) then
setenv G4VIS_BUILD_OIX_DRIVER 1
echo "On this machine the G4VIS_BUILD_OIX_DRIVER=\$G4VIS_BUILD_OIX_DRIVER"
endif

#+
if ( X$g4vis_build_raytracerx_driver == Xy ) then
setenv G4VIS_BUILD_RAYTRACERX_DRIVER 1
echo "On this machine the G4VIS_BUILD_RAYTRACERX_DRIVER=\$G4VIS_BUILD_RAYTRACERX_DRIVER"
endif

#+
if ( X$g4vis_build_oiwin32_driver == Xy ) then
setenv G4VIS_BUILD_OIWIN32_DRIVER 1
echo "On this machine the G4VIS_BUILD_OIWIN32_DRIVER=\$G4VIS_BUILD_OIWIN32_DRIVER"
endif

#+
if ( X$g4vis_build_vrml_driver == Xy ) then
setenv G4VIS_BUILD_VRML_DRIVER 1
echo "On this machine the G4VIS_BUILD_VRML_DRIVER=\$G4VIS_BUILD_VRML_DRIVER"
endif

#+
if ( X$g4vis_use_dawn == Xy ) then
setenv G4VIS_USE_DAWN 1
echo "On this machine the G4VIS_USE_DAWN=\$G4VIS_USE_DAWN"
endif

#+
if ( X$g4vis_use_openglx == Xy ) then
setenv G4VIS_USE_OPENGLX 1
echo "On this machine the G4VIS_USE_OPENGLX=\$G4VIS_USE_OPENGLX"
endif

#+
if ( X$g4vis_use_openglxm == Xy ) then
setenv G4VIS_USE_OPENGLXM 1
echo "On this machine the G4VIS_USE_OPENGLXM=\$G4VIS_USE_OPENGLXM"
endif

#+
if ( X$g4vis_use_openglwin32 == Xy ) then
setenv G4VIS_USE_OPENGLWIN32 1
echo "On this machine the G4VIS_USE_OPENGLWIN32=\$G4VIS_USE_OPENGLWIN32"
endif

#+
if ( X$g4vis_use_oix == Xy ) then
setenv G4VIS_USE_OIX 1
echo "On this machine the G4VIS_USE_OIX=\$G4VIS_USE_OIX"
endif

#+
if ( X$g4vis_use_raytracerx == Xy ) then
setenv G4VIS_USE_RAYTRACERX 1
echo "On this machine the G4VIS_USE_RAYTRACERX=\$G4VIS_USE_RAYTRACERX"
endif

#+
if ( X$g4vis_use_oiwin32 == Xy ) then
setenv G4VIS_USE_OIWIN32 1
echo "On this machine the G4VIS_USE_OIWIN32=\$G4VIS_USE_OIWIN32"
endif

#+
if ( X$g4vis_use_vrml == Xy ) then
setenv G4VIS_USE_VRML 1
echo "On this machine the G4VIS_USE_VRML=\$G4VIS_USE_VRML"
endif

#+
if ( X$g4vis_oglhome != X )  then
setenv OGLHOME "$g4vis_oglhome"
echo "On this machine the OGLHOME=\$OGLHOME"
endif 

#+
if ( X$g4vis_oivhome != X )  then
setenv OIVHOME "$g4vis_oivhome"
echo "On this machine the OIVHOME=\$OIVHOME"
endif 

#+
if ( X$g4vis_build_openglxm_driver != X )  then
setenv XMFLAGS "$g4vis_xmflags"
echo "On this machine the XMFLAGS=\$XMFLAGS"
endif 

#+
if ( X$g4vis_build_openglxm_driver != X )  then
setenv XMLIBS "$g4vis_xmlibs"
echo "On this machine the XMLIBS=\$XMLIBS"
endif 

#+
if ( X$g4ui_build_xm_session != X )  then
setenv XMFLAGS "$g4vis_xmflags"
echo "On this machine the XMFLAGS=\$XMFLAGS"
endif 

#+
if ( X$g4ui_build_xm_session != X )  then
setenv XMLIBS "$g4vis_xmlibs"
echo "On this machine the XMLIBS=\$XMLIBS"
endif 

#+
if ( X$g4ui_build_xaw_session != X )  then
setenv XAWFLAGS "$g4ui_xawflags"
echo "On this machine the XAWFLAGS=\$XAWFLAGS"
endif 

#+
if ( X$g4ui_build_xaw_session != X )  then
setenv XAWLIBS "$g4ui_xawlibs"
echo "On this machine the XAWLIBS=\$XAWLIBS"
endif 

#
# Use G3TOG4 module
#
#+
if ( X$g4wlib_build_g3tog4 == Xy )  then
setenv G4LIB_BUILD_G3TOG4 1
echo "On this machine the G4LIB_BUILD_G3TOG4=\$G4LIB_BUILD_G3TOG4"
endif 

if ( X$g4w_use_g3tog4 == Xy )  then
setenv G4LIB_USE_G3TOG4 1
echo "On this machine the G4LIB_USE_G3TOG4=\$G4LIB_USE_G3TOG4"
endif 


#
# Use ZLIB module
#
#+
if ( X$g4wlib_build_zlib == Xy )  then
setenv G4LIB_BUILD_ZLIB 1
echo "On this machine the G4LIB_BUILD_ZLIB=\$G4LIB_BUILD_ZLIB"
endif 

if ( X$g4w_use_ZLIB == Xy )  then
setenv G4LIB_USE_ZLIB 1
echo "On this machine the G4LIB_USE_ZLIB=\$G4LIB_USE_ZLIB"
endif 

#+
#
# g4shared
#
if ( X$g4lib_build_shared == Xy ) then
setenv G4LIB_BUILD_SHARED 1
echo "On this machine the G4LIB_BUILD_SHARED=\$G4LIB_BUILD_SHARED"
endif

if ( X$g4lib_build_static == Xy ) then
setenv G4LIB_BUILD_STATIC 1
echo "On this machine the G4LIB_BUILD_STATIC=\$G4LIB_BUILD_STATIC"
endif

#+
#
# g4granular
#
if ( X$g4lib_use_granular == Xy ) then
setenv G4LIB_USE_GRANULAR 1
echo "On this machine the G4LIB_USE_GRANULAR=\$G4LIB_USE_GRANULAR"
endif

# Check for Windows!
if ( "X\$G4SYSTEM" != "XWIN32-VC" && "X\$G4SYSTEM" != "XWIN32-VC7" ) then
setenv G4UI_USE_TCSH 1
echo "On this machine the G4UI_USE_TCSH=\$G4UI_USE_TCSH"
endif

#####################################################################

!GROK!THIS!

: In the following dollars and backticks do not need the extra backslash.
$spitshell >>env.csh <<'!NO!SUBS!'


#+
#
# G4WORKDIR
#
if ( ${?G4WORKDIR} ) then
echo "In your environment you have the G4WORKDIR=$G4WORKDIR"
else
# Check for Windows!
if ( "X$G4SYSTEM" == "XWIN32-VC" || "X$G4SYSTEM" == "XWIN32-VC7" ) then

echo "G4WORKDIR will be set to c:/geant4 (in "native" Windows syntax)."
setenv G4WORKDIR "c:/geant4"

else # if Unix

echo "G4WORKDIR will be set to HOME/geant4=$HOME/geant4"
setenv G4WORKDIR $HOME/geant4

endif
endif

#
# Shared libraries
#
if ( ${?G4LIB_BUILD_SHARED} ) then

#
# Warning about LD_LIBRARY_PATH
#
echo ''
echo 'LD_LIBRARY_PATH is set to include CLHEP and Geant4 libraries. '
echo ''
echo 'NOTE : verify that the correct path for the CLHEP library'
echo '- $CLHEP_BASE_DIR/lib - is included in the $LD_LIBRARY_PATH'
echo 'variable in your environment, and no other outdated installations'
echo 'of CLHEP are by chance referenced through $LD_LIBRARY_PATH and'
echo 'present in your system !'
echo ''

if ( ${?LD_LIBRARY_PATH} )  then
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${G4LIB}/${G4SYSTEM}
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${CLHEP_BASE_DIR}/lib
else
setenv LD_LIBRARY_PATH ${G4LIB}/${G4SYSTEM}
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${CLHEP_BASE_DIR}/lib
endif
endif

#
# Shared libraries for Mac
#
if ( "X$G4SYSTEM" == "XDarwin-g++" ) then

if ( ${?G4LIB_BUILD_SHARED} ) then
if ( ${?DYLD_LIBRARY_PATH} )  then
setenv DYLD_LIBRARY_PATH ${DYLD_LIBRARY_PATH}:${G4LIB}/${G4SYSTEM}
setenv DYLD_LIBRARY_PATH ${DYLD_LIBRARY_PATH}:${CLHEP_BASE_DIR}/lib
else
setenv DYLD_LIBRARY_PATH ${G4LIB}/${G4SYSTEM}
setenv DYLD_LIBRARY_PATH ${DYLD_LIBRARY_PATH}:${CLHEP_BASE_DIR}/lib
endif
endif

endif

#
# For BIN
#
setenv PATH ${PATH}:${G4WORKDIR}/bin/${G4SYSTEM}

!NO!SUBS!
chmod 755 env.csh
$eunicefix env.csh
