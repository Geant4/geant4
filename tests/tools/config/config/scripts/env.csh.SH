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
# g4system.U
#
#+
setenv G4SYSTEM "$g4system"
echo "On this machine the G4SYSTEM=\$G4SYSTEM"

#
# g4dirs.U
#
#+
setenv G4INSTALL "$g4install"
echo "On this machine the G4INSTALL=\$G4INSTALL"

#+
#?#G4INCLUDE="$g4include"

#+
if ( X$g4lib != X ) then
setenv G4LIB "$g4lib"
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
setenv NeutronHPCrossSections "$g4neutronhpcrosssections"
echo "On this machine the NeutronHPCrossSections=\$NeutronHPCrossSections"
endif

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
if ( X$g4debug != X ) then
setenv G4DEBUG "$g4debug"
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

#+
#
# g4persistency
#
if ( X$g4use_hepodbms == Xy ) then
setenv G4USE_HEPODBMS 1
source  \$G4INSTALL/examples/extended/persistency/PersistentEx01/g4odbms_setup.csh
echo "On this machine the G4USE_HEPODBMS=\$G4USE_HEPODBMS"
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

#+
if ( X$g4ui_build_wo_session == Xy ) then
setenv G4UI_BUILD_WO_SESSION 1
echo "On this machine the G4UI_BUILD_WO_SESSION=\$G4UI_BUILD_WO_SESSION"
endif

#+
if ( X$g4ui_use_wo == Xy ) then
setenv G4UI_USE_WO 1
echo "On this machine the G4UI_USE_WO=\$G4UI_USE_WO"
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
if ( X$g4vis_build_opacs_driver == Xy ) then
setenv G4VIS_BUILD_OPACS_DRIVER 1
echo "On this machine the G4VIS_BUILD_OPACS_DRIVER=\$G4VIS_BUILD_OPACS_DRIVER"
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
if ( X$g4vis_use_opacs == Xy ) then
setenv G4VIS_USE_OPACS 1
echo "On this machine the G4VIS_USE_OPACS=\$G4VIS_USE_OPACS"
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

if ( X$g4lib_use_shared == Xy ) then
setenv G4LIB_USE_SHARED 1
echo "On this machine the G4LIB_USE_SHARED=\$G4LIB_USE_SHARED"
endif

if ( X$g4lib_use_static == Xy ) then
setenv G4LIB_USE_STATIC 1
echo "On this machine the G4LIB_USE_STATIC=\$G4LIB_USE_STATIC"
endif

#+
#
# g4granular
#
if ( X$g4lib_use_granular == Xy ) then
setenv G4LIB_USE_GRANULAR 1
echo "On this machine the G4LIB_USE_GRANULAR=\$G4LIB_USE_GRANULAR"
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
echo "G4WORKDIR will be set to HOME=$HOME."
setenv G4WORKDIR $HOME
endif


!NO!SUBS!
chmod 755 env.csh
$eunicefix env.csh
