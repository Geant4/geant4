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
#echo "Extracting install.sh (with variable substitutions)"
: This section of the file will have variable substitutions done on it.
: Move anything that needs config subs from !NO!SUBS! section to !GROK!THIS!.
: Protect any dollar signs and backticks that you do not want interpreted
: by putting a backslash in front.  You may delete these comments.
$spitshell >install.sh <<!GROK!THIS!
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

if ( X$g4base != X ) then
setenv G4BASE "$g4base"
echo "On this machine the G4BASE=\$G4BASE"
endif

if ( X$g4workdir != X ) then
setenv G4WORKDIR "$g4workdir"
echo "On this machine the G4WORKDIR=\$G4WORKDIR"
endif

#+
if ( X$g4tmp != X ) then
setenv G4TMP "$g4tmp"
echo "On this machine the G4TMP=\$G4TMP"
endif
#+
if ( X$g4lib != X ) then
setenv G4LIB "$g4lib"
echo "On this machine the G4LIB=\$G4LIB"
endif
#+
if ( X$g4bin != X ) then
setenv G4BIN "$g4bin"
echo "On this machine the G4BIN=\$G4BIN"
endif

#+
if ( X$g4data != X ) then
setenv G4DATA "$g4data"
echo "On this machine the G4DATA=\$G4DATA"
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
echo "On this machine the G4LEVELGAMMADATA=\$G4LEVELGAMMADATA"
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
.  \$G4INSTALL/examples/extended/persistency/PersistentEx01/g4odbms_setup.sh
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

#+ set by default in U file
if ( X$g4ui_build_terminal_session == Xy ) then
setenv G4UI_BUILD_TERMINAL_SESSION 1
echo "On this machine the G4UI_BUILD_TERMINAL_SESSION=\$G4UI_BUILD_TERMINAL_SESSION"
endif

#+ set by default in U file
if ( X$g4ui_use_terminal == Xy ) then
setenv G4UI_USE_TERMINAL 1
echo "On this machine the G4UI_USE_TERMINAL=\$G4UI_USE_TERMINAL"
endif

#+ set by default in U file
if ( X$g4ui_build_gag_session == Xy ) then
setenv G4UI_BUILD_GAG_SESSION 1
echo "On this machine the G4UI_BUILD_GAG_SESSION=\$G4UI_BUILD_GAG_SESSION"
endif

#+ set by default in U file
if ( X$g4ui_use_gag == Xy ) then
setenv G4UI_USE_GAG 1
echo "On this machine the G4UI_USE_GAG=\$G4UI_USE_GAG"
endif

#+ set by default in U file
if ( X$g4ui_use_tcsh == Xy ) then
setenv G4UI_USE_TCSH 1
echo "On this machine the G4UI_USE_TCSH=\$G4UI_USE_TCSH"
endif

#
# g4vis
#
#+
if ( X$g4vis_build_dawn_driver == Xy ) then
setenv G4VIS_BUILD_DAWN_DRIVER 1
echo "On this machine the G4VIS_BUILD_DAWN_DRIVER=\$G4VIS_BUILD_DAWN_DRIVER"
endif

#+ set by default in U file
if ( X$g4vis_build_dawnfile_driver == Xy ) then
setenv G4VIS_BUILD_DAWNFILE_DRIVER 1
echo "On this machine the G4VIS_BUILD_DAWNFILE_DRIVER=\$G4VIS_BUILD_DAWNFILE_DRIVER"
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
if ( X$g4vis_build_oix_driver == Xy ) then
setenv G4VIS_BUILD_OIX_DRIVER 1
echo "On this machine the G4VIS_BUILD_OIX_DRIVER=\$G4VIS_BUILD_OIX_DRIVER"
endif

#+ set by default in U file
if ( X$g4vis_build_raytracer_driver == Xy ) then
setenv G4VIS_BUILD_RAYTRACER_DRIVER 1
echo "On this machine the G4VIS_BUILD_RAYTRACER_DRIVER=\$G4VIS_BUILD_RAYTRACER_DRIVER"
endif

#+
if ( X$g4vis_build_vrml_driver == Xy ) then
setenv G4VIS_BUILD_VRML_DRIVER 1
echo "On this machine the G4VIS_BUILD_VRML_DRIVER=\$G4VIS_BUILD_VRML_DRIVER"
endif

#+ set by default in U file
if ( X$g4vis_build_vrmlfile_driver == Xy ) then
setenv G4VIS_BUILD_VRMLFILE_DRIVER 1
echo "On this machine the G4VIS_BUILD_VRMLFILE_DRIVER=\$G4VIS_BUILD_VRMLFILE_DRIVER"
endif

# set by default in U file
if ( X$g4vis_build_asciitree_driver == Xy ) then
setenv G4VIS_BUILD_ASCIITREE_DRIVER 1
echo "On this machine the G4VIS_BUILD_ASCIITREE_DRIVER=\$G4VIS_BUILD_ASCIITREE_DRIVER"
endif

#+
if ( X$g4vis_use_dawn == Xy ) then
setenv G4VIS_USE_DAWN 1
echo "On this machine the G4VIS_USE_DAWN=\$G4VIS_USE_DAWN"
endif

#+ set by default in U file
if ( X$g4vis_use_dawnfile == Xy ) then
setenv G4VIS_USE_DAWNFILE 1
echo "On this machine the G4VIS_USE_DAWNFILE=\$G4VIS_USE_DAWNFILE"
endif

#+
if ( X$g4vis_use_opacs == Xy ) then
setenv G4VIS_USE_OPACS 1
echo "On this machine the G4VIS_USE_OPACS=\$G4VIS_USE_OPACS"
endif

#?
if ( X$g4vis_use_opengl == Xy ) then
setenv G4VIS_USE_OPENGL 1
else
setenv G4VIS_USE_OPENGL ""
endif

#+
if ( X$g4vis_use_OPENGLX == Xy ) then
setenv G4VIS_USE_OPENGLX 1
echo "On this machine the G4VIS_USE_OPENGLX=\$G4VIS_USE_OPENGLX"
endif

#+
if ( X$g4vis_use_oix == Xy ) then
setenv G4VIS_USE_OIX 1
echo "On this machine the G4VIS_USE_OIX=\$G4VIS_USE_OIX"
endif

#+ set by default in U file
if ( X$g4vis_use_raytracer == Xy ) then
setenv G4VIS_USE_RAYTRACER 1
echo "On this machine the G4VIS_USE_RAYTRACER=\$G4VIS_USE_RAYTRACER"
endif

#+
if ( X$g4vis_use_vrml == Xy ) then
setenv G4VIS_USE_VRML 1
echo "On this machine the G4VIS_USE_VRML=\$G4VIS_USE_VRML"
endif

#+ set by default in U file
if ( X$g4vis_use_vrmlfile == Xy ) then
setenv G4VIS_USE_VRMLFILE 1
echo "On this machine the G4VIS_USE_VRMLFILE=\$G4VIS_USE_VRMLFILE"
endif

#+ set by default in U file
if ( X$g4vis_use_asciitree == Xy ) then
setenv G4VIS_USE_ASCIITREE 1
echo "On this machine the G4VIS_USE_ASCIITREE=\$G4VIS_USE_ASCIITREE"
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
#export G4SYSTEM

#export G4INSTALL
#?#export G4INCLUDE
#export G4BASE
#export G4WORKDIR

#export G4TMP
#export G4LIB
#export G4BIN

#export G4DATA
#export G4LEVELGAMMADATA

#export CLHEP_BASE_DIR
#export CLHEP_INCLUDE_DIR
#export CLHEP_LIB_DIR
#export CLHEP_LIB

#export G4DEBUG

#export G4ANALYSIS_USE

#export G4USE_HEPODBMS

#export G4UI_BUILD_XAW_SESSION
#export G4UI_USE_XAW
#export G4UI_BUILD_XM_SESSION
#export G4UI_USE_XM
#export G4UI_BUILD_WIN32_SESSION
#export G4UI_USE_WIN32

#export G4UI_BUILD_TERMINAL_SESSION
#export G4UI_USE_TERMINAL

#export G4UI_BUILD_GAG_SESSION
#export G4UI_USE_GAG

#export G4UI_USE_TCSH

#export G4VIS_BUILD_DAWN_DRIVER
#export G4VIS_USE_DAWN
#export G4VIS_BUILD_DAWNFILE_DRIVER
#export G4VIS_USE_DAWNFILE
#export G4VIS_BUILD_OPACS_DRIVER
#export G4VIS_USE_OPACS
#export G4VIS_BUILD_OPENGLX_DRIVER
#export G4VIS_USE_OPENGLX
#export G4VIS_BUILD_OPENGLXM_DRIVER
#export G4VIS_USE_OPENGLXM
#export G4VIS_BUILD_OIX_DRIVER
#export G4VIS_USE_OIX
#export G4VIS_BUILD_RAYTRACER_DRIVER
#export G4VIS_USE_RAYTRACER
#export G4VIS_BUILD_VRML_DRIVER
#export G4VIS_USE_VRML
#export G4VIS_BUILD_VRMLFILE_DRIVER
#export G4VIS_USE_VRMLFILE
#export G4VIS_BUILD_ASCIITREE_DRIVER
#export G4VIS_USE_ASCIITREE

#export G4LIB_BUILD_SHARED
#export G4LIB_BUILD_STATIC
#export G4LIB_USE_SHARED
#export G4LIB_USE_STATIC


###############################################


#####################################################################
#echo "On this machine the G4SYSTEM=\$G4SYSTEM"

#echo "On this machine the G4INSTALL=\$G4INSTALL"
#?#echo "On this machine the G4INCLUDE=\$G4INCLUDE"
#echo "On this machine the G4BASE=\$G4BASE"
#echo "On this machine the G4WORKDIR=\$G4WORKDIR"
#echo "On this machine the G4TMP=\$G4TMP"
#echo "On this machine the G4LIB=\$G4LIB"
#echo "On this machine the G4BIN=\$G4BIN"

#echo "On this machine the G4DATA=\$G4DATA"
#echo "On this machine the G4LEVELGAMMADATA=\$G4LEVELGAMMADATA"

#echo "On this machine the CLHEP_BASE_DIR=\$CLHEP_BASE_DIR"
#echo "On this machine the CLHEP_INCLUDE_DIR=\$CLHEP_INCLUDE_DIR"
#echo "On this machine the CLHEP_LIB_DIR=\$CLHEP_LIB_DIR"
#echo "On this machine the CLHEP_LIB=\$CLHEP_LIB"

#echo "On this machine the G4DEBUG=\$G4DEBUG"

#echo "On this machine the G4UI_BUILD_TERMINAL_SESSION=\$G4UI_BUILD_TERMINAL_SESSION"
#echo "On this machine the G4UI_USE_TERMINAL=\$G4UI_USE_TERMINAL"
#echo "On this machine the G4UI_BUILD_GAG_SESSION=\$G4UI_BUILD_GAG_SESSION"
#echo "On this machine the G4UI_USE_GAG=\$G4UI_USE_GAG"
#echo "On this machine the G4UI_BUILD_XAW_SESSION=\$G4UI_BUILD_XAW_SESSION"
#echo "On this machine the G4UI_USE_XAW=\$G4UI_USE_XAW"
#echo "On this machine the G4UI_BUILD_XM_SESSION=\$G4UI_BUILD_XM_SESSION"
#echo "On this machine the G4UI_USE_XM=\$G4UI_USE_XM"
#echo "On this machine the G4UI_BUILD_WIN32_SESSION=\$G4UI_BUILD_WIN32_SESSION"
#echo "On this machine the G4UI_USE_WIN32=\$G4UI_USE_WIN32"
#echo "On this machine the G4UI_USE_TCSH=\$G4UI_USE_TCSH"

#echo "On this machine the G4ANALYSIS_USE=\$G4ANALYSIS_USE"

#echo "On this machine the G4USE_HEPODBMS=\$G4USE_HEPODBMS"

#echo "On this machine the G4VIS_BUILD_DAWN_DRIVER=\$G4VIS_BUILD_DAWN_DRIVER"
#echo "On this machine the G4VIS_USE_DAWN=\$G4VIS_USE_DAWN"
#echo "On this machine the G4VIS_BUILD_DAWNFILE_DRIVER=\$G4VIS_BUILD_DAWNFILE_DRIVER"
#echo "On this machine the G4VIS_USE_DAWNFILE=\$G4VIS_USE_DAWNFILE"
#echo "On this machine the G4VIS_BUILD_OPACS_DRIVER=\$G4VIS_BUILD_OPACS_DRIVER"
#echo "On this machine the G4VIS_USE_OPACS=\$G4VIS_USE_OPACS"
#echo "On this machine the G4VIS_BUILD_OPENGLX_DRIVER=\$G4VIS_BUILD_OPENGLX_DRIVER"

echo "On this machine the G4VIS_USE_OPENGLX=\$G4VIS_USE_OPENGLX"

#echo "On this machine the G4VIS_BUILD_OPENGLXM_DRIVER=\$G4VIS_BUILD_OPENGLXM_DRIVER"
#echo "On this machine the G4VIS_USE_OPENGLXM=\$G4VIS_USE_OPENGLXM"
#echo "On this machine the G4VIS_BUILD_OIX_DRIVER=\$G4VIS_BUILD_OIX_DRIVER"
#echo "On this machine the G4VIS_USE_OIX=\$G4VIS_USE_OIX"
#echo "On this machine the G4VIS_BUILD_RAYTRACER_DRIVER=\$G4VIS_BUILD_RAYTRACER_DRIVER"
#echo "On this machine the G4VIS_USE_RAYTRACER=\$G4VIS_USE_RAYTRACER"
#echo "On this machine the G4VIS_BUILD_VRML_DRIVER=\$G4VIS_BUILD_VRML_DRIVER"
#echo "On this machine the G4VIS_USE_VRML=\$G4VIS_USE_VRML"
#echo "On this machine the G4VIS_BUILD_VRMLFILE_DRIVER=\$G4VIS_BUILD_VRMLFILE_DRIVER"
#echo "On this machine the G4VIS_USE_VRMLFILE=\$G4VIS_USE_VRMLFILE"
#echo "On this machine the G4VIS_BUILD_ASCIITREE_DRIVER=\$G4VIS_BUILD_ASCIITREE_DRIVER"
#echo "On this machine the G4VIS_USE_ASCIITREE=\$G4VIS_USE_ASCIITREE"

#echo "On this machine the G4LIB_BUILD_SHARED=\$G4LIB_BUILD_SHARED"
#echo "On this machine the G4LIB_BUILD_STATIC=\$G4LIB_BUILD_STATIC"
#echo "On this machine the G4LIB_USE_SHARED=\$G4LIB_USE_SHARED"
#echo "On this machine the G4LIB_USE_STATIC=\$G4LIB_USE_STATIC"

#echo "On this machine the G4LIB_USE_GRANULAR=\$G4LIB_USE_GRANULAR"

###############################################################

!GROK!THIS!

: In the following dollars and backticks do not need the extra backslash.
$spitshell >>install.sh <<'!NO!SUBS!'
!NO!SUBS!
chmod 755 install.sh
$eunicefix install.sh
