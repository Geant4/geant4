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
set G4SYSTEM="$g4system"

#
# g4dirs.U
#
set G4INSTALL="$g4install"
set G4INCLUDE="$g4include"
set G4BASE="$g4base"
set G4WORKDIR="$g4workdir"
set G4TMP="$g4tmp"
set G4LIB="$g4lib"
set G4BIN="$g4bin"
set G4DATA="$g4data"
set G4LEVELGAMMADATA="$g4levelgammadata"
#+
set G4RADIOACTIVEDATA="$g4radioactivedata"
#+
set G4LEDATA="$g4ledata"
#+
set NeutronHPCrossSections="$g4neutronhpcrosssections"

#
# g4clhep.U
#
set CLHEP_BASE_DIR="$g4clhep_base_dir"
set CLHEP_INCLUDE_DIR="$g4clhep_include_dir"
set set CLHEP_LIB_DIR="$g4clhep_lib_dir"
set CLHEP_LIB="$g4clhep_lib"

#
# g4ospace
#
if ( X$g4use_ospace == Xy ) then
set G4USE_OSPACE=1
else
set G4USE_OSPACE=""
endif

set G4OSPACE_BASE_DIR="$g4ospace_base_dir"

#
# g4debug
#
set G4DEBUG="$g4debug"

#
# g4analysis
#
if ( X$g4wanalysis_build == Xy ) then
set G4ANALYSIS_BUILD=1
else
set G4ANALYSIS_BUILD=""
endif 

if ( X$g4wanalysis_build_lab == Xy ) then
set G4ANALYSIS_BUILD_LAB=1
else
set G4ANALYSIS_BUILD_LAB=""
endif 

if ( X$g4wanalysis_build_jas == Xy ) then
set G4ANALYSIS_BUILD_JAS=1
else
set G4ANALYSIS_BUILD_JAS=""
endif 

if ( X$g4wanalysis_build_lizard == Xy ) then
set G4ANALYSIS_BUILD_LIZARD=1
else
set G4ANALYSIS_BUILD_LIZARD=""
endif 

if ( X$g4wanalysis_use == Xy ) then
set G4ANALYSIS_USE=1
else
set G4ANALYSIS_USE=""
endif 

if ( X$g4wanalysis_use_lab == Xy ) then
set G4ANALYSIS_USE_LAB=1
else
set G4ANALYSIS_USE_LAB=""
endif 

if ( X$g4wanalysis_use_jas == Xy ) then
set G4ANALYSIS_USE_JAS=1
else
set G4ANALYSIS_USE_JAS=""
endif 

if ( X$g4wanalysis_use_lizard == Xy ) then
set G4ANALYSIS_USE_LIZARD=1
else
set G4ANALYSIS_USE_LIZARD=""
endif 

#+
#
# g4persistency
#
if ( X$g4use_hepodbms == Xy ) then
set G4USE_HEPODBMS=1
source  \$G4INSTALL/examples/extended/persistency/PersistentEx01/g4odbms_setup.csh
else
set G4USE_HEPODBMS=""
endif 

#
# g4ui
#
if ( X$g4ui_build_terminal_session == Xy ) then
set G4UI_BUILD_TERMINAL_SESSION=1
else
set G4UI_BUILD_TERMINAL_SESSION=""
endif 

if ( X$g4ui_use_terminal == Xy ) then
set G4UI_USE_TERMINAL=1
else
set G4UI_USE_TERMINAL=""
endif 

if ( X$g4ui_build_gag_session == Xy ) then
set G4UI_BUILD_GAG_SESSION=1
else
set G4UI_BUILD_GAG_SESSION=""
endif 

if ( X$g4ui_use_gag == Xy ) then
set G4UI_USE_GAG=1
else
set G4UI_USE_GAG=""
endif 

if ( X$g4ui_build_xaw_session == Xy ) then
set G4UI_BUILD_XAW_SESSION=1
else
set G4UI_BUILD_XAW_SESSION=""
endif 

if ( X$g4ui_use_xaw == Xy ) then
set G4UI_USE_XAW=1
else
set G4UI_USE_XAW=""
endif 

if ( X$g4ui_build_xm_session == Xy ) then
set G4UI_BUILD_XM_SESSION=1
else
set G4UI_BUILD_XM_SESSION=""
endif 

if ( X$g4ui_use_xm == Xy ) then
set G4UI_USE_XM=1
else
set G4UI_USE_XM=""
endif 

if ( X$g4ui_build_win32_session == Xy ) then
set G4UI_BUILD_WIN32_SESSION=1
else
set G4UI_BUILD_WIN32_session=""
endif 

if ( X$g4ui_use_win32 == Xy ) then
set G4UI_USE_WIN32=1
else
set G4UI_USE_WIN32=""
endif 

if ( X$g4ui_use_tcsh == Xy ) then
set G4UI_USE_TCSH=1
else
set G4UI_USE_TCSH=""
endif 

#
# g4vis
#
if ( X$g4vis_build_dawn_driver == Xy ) then
set G4VIS_BUILD_DAWN_DRIVER=1
else
set G4VIS_BUILD_DAWN_DRIVER=""
endif 

if ( X$g4vis_build_dawnfile_driver == Xy ) then
set G4VIS_BUILD_DAWNFILE_DRIVER=1
else
set G4VIS_BUILD_DAWN_DRIVER=""
endif 

if ( X$g4vis_build_opacs_driver == Xy ) then
set G4VIS_BUILD_OPACS_DRIVER=1
else
set G4VIS_BUILD_OPACS_DRIVER=""
endif 

if ( X$g4vis_build_openglx_driver == Xy ) then
set G4VIS_BUILD_OPENGLX_DRIVER=1
else
set G4VIS_BUILD_OPENGLX_DRIVER=""
endif 

if ( X$g4vis_build_openglxm_driver == Xy ) then
set G4VIS_BUILD_OPENGLXM_DRIVER=1
else
set G4VIS_BUILD_OPENGLXM_DRIVER=""
endif 

if ( X$g4vis_build_oix_driver == Xy ) then
set G4VIS_BUILD_OIX_DRIVER=1
else
set G4VIS_BUILD_OIX_DRIVER=""
endif 

if ( X$g4vis_build_raytracer_driver == Xy ) then
set G4VIS_BUILD_RAYTRACER_DRIVER=1
else
set G4VIS_BUILD_RAYTRACER_DRIVER=""
endif 

if ( X$g4vis_build_vrml_driver == Xy ) then
set G4VIS_BUILD_VRML_DRIVER=1
else
set G4VIS_BUILD_VRML_DRIVER=""
endif 

if ( X$g4vis_build_vrmlfile_driver == Xy ) then
set G4VIS_BUILD_VRMLFILE_DRIVER=1
else
set G4VIS_BUILD_VRMLFILE_DRIVER=""
endif 

if ( X$g4vis_build_asciitree_driver == Xy ) then
set G4VIS_BUILD_ASCIITREE_DRIVER=1
else
set G4VIS_BUILD_ASCIITREE_DRIVER=""
endif 

if ( X$g4vis_use_dawn == Xy ) then
set G4VIS_USE_DAWN=1
else
set G4VIS_USE_DAWN=""
endif 

if ( X$g4vis_use_dawnfile == Xy ) then
set G4VIS_USE_DAWNFILE=1
else
set G4VIS_USE_DAWN=""
endif 

if ( X$g4vis_use_opacs == Xy ) then
set G4VIS_USE_OPACS=1
else
set G4VIS_USE_OPACS=""
endif

if ( X$g4vis_use_opengl == Xy ) then
set G4VIS_USE_OPENGL=1
else
set G4VIS_USE_OPENGL=""
endif 

if ( X$g4vis_use_OPENGLX == Xy ) then
set G4VIS_USE_OPENGLX=1
else
set G4VIS_USE_OPENGLX=""
endif 

if ( X$g4vis_use_oix == Xy ) then
set G4VIS_USE_OIX=1
else
set G4VIS_USE_OIX=""
endif 

if ( X$g4vis_use_raytracer == Xy ) then
set G4VIS_USE_RAYTRACER=1
else
set G4VIS_USE_RAYTRACER=""
endif 

if ( X$g4vis_use_vrml == Xy ) then
set G4VIS_USE_VRML=1
else
set G4VIS_USE_VRML=""
endif 

if ( X$g4vis_use_vrmlfile == Xy ) then
set G4VIS_USE_VRMLFILE=1
else
set G4VIS_USE_VRMLFILE=""
endif 

if ( X$g4vis_use_asciitree == Xy ) then
set G4VIS_USE_ASCIITREE=1
else
set G4VIS_USE_ASCIITREE=""
endif

#
# g4shared
#
set G4LIB_BUILD_SHARED="$g4lib_build_shared"
set G4LIB_BUILD_STATIC="$g4lib_build_static"
set G4LIB_USE_SHARED="$g4lib_use_shared"
set G4LIB_USE_STATIC="$g4lib_use_static"

if ( X$g4lib_build_shared == Xy ) then
set G4LIB_BUILD_SHARED=1
else
set G4LIB_BUILD_SHARED=""
endif 

if ( X$g4lib_build_static == Xy ) then
set G4LIB_BUILD_STATIC=1
else
set G4LIB_BUILD_STATIC=""
endif 

if ( X$g4lib_use_shared == Xy ) then
set G4LIB_USE_SHARED=1
else
set G4LIB_USE_SHARED=""
endif 

if ( X$g4lib_use_static == Xy ) then
set G4LIB_USE_STATIC=1
else
set G4LIB_USE_STATIC=""
endif 

#
# g4granular
#
if ( X$g4lib_use_granular == Xy ) then
set G4LIB_USE_GRANULAR=1
else
set G4LIB_USE_GRANULAR=""
endif 

#####################################################################
setenv G4SYSTEM	 $G4SYSTEM

setenv G4INSTALL	$G4INSTALL
setenv G4INCLUDE	 $G4INCLUDE
setenv G4BASE	 $G4BASE
setenv G4WORKDIR	 $G4WORKDIR
setenv G4TMP	 $G4TMP
setenv G4LIB	 $G4LIB
setenv G4BIN	 $G4BIN
setenv G4DATA	 $G4DATA
setenv G4LEVELGAMMADATA	 $G4LEVELGAMMADATA

setenv CLHEP_BASE_DIR	 $CLHEP_BASE_DIR
setenv CLHEP_INCLUDE_DIR	 $CLHEP_INCLUDE_DIR
setenv CLHEP_LIB_DIR	 $CLHEP_LIB_DIR
setenv CLHEP_LIB	 $CLHEP_LIB

setenv G4USE_OSPACE	 $G4USE_OSPACE
setenv G4OSPACE_BASE_DIR	 $G4OSPACE_BASE_DIR

setenv G4LIB_USE_GRANULAR	 $G4LIB_USE_GRANULAR

setenv G4DEBUG	 $G4DEBUG

setenv G4ANALYSIS_BUILD	 $G4ANALYSIS_BUILD
setenv G4ANALYSIS_BUILD_LAB	 $G4ANALYSIS_BUILD_LAB
setenv G4ANALYSIS_BUILD_JAS	 $G4ANALYSIS_BUILD_JAS
setenv G4ANALYSIS_BUILD_LIZARD	 $G4ANALYSIS_BUILD_LIZARD
setenv G4ANALYSIS_USE	 $G4ANALYSIS_USE
setenv G4ANALYSIS_USE_LAB	 $G4ANALYSIS_USE_LAB
setenv G4ANALYSIS_USE_JAS	 $G4ANALYSIS_USE_JAS
setenv G4ANALYSIS_USE_LIZARD	 $G4ANALYSIS_USE_LIZARD

setenv G4USE_HEPODBMS	 $G4USE_HEPODBMS

setenv G4UI_BUILD_TERMINAL_SESSION	 $G4UI_BUILD_TERMINAL_SESSION
setenv G4UI_USE_TERMINAL	 $G4UI_USE_TERMINAL
setenv G4UI_BUILD_GAG_SESSION	 $G4UI_BUILD_GAG_SESSION
setenv G4UI_USE_GAG	 $G4UI_USE_GAG
setenv G4UI_BUILD_XAW_SESSION	 $G4UI_BUILD_XAW_SESSION
setenv G4UI_USE_XAW	 $G4UI_USE_XAW
setenv G4UI_BUILD_XM_SESSION	 $G4UI_BUILD_XM_SESSION
setenv G4UI_USE_XM	 $G4UI_USE_XM
setenv G4UI_BUILD_WIN32_SESSION	 $G4UI_BUILD_WIN32_SESSION
setenv G4UI_USE_WIN32	 $G4UI_USE_WIN32
setenv G4UI_USE_TCSH	 $G4UI_USE_TCSH

setenv G4VIS_BUILD_DAWN_DRIVER	 $G4VIS_BUILD_DAWN_DRIVER
setenv G4VIS_USE_DAWN	 $G4VIS_USE_DAWN
setenv G4VIS_BUILD_DAWNFILE_DRIVER	 $G4VIS_BUILD_DAWNFILE_DRIVER
setenv G4VIS_USE_DAWNFILE	 $G4VIS_USE_DAWNFILE
setenv G4VIS_BUILD_OPACS_DRIVER	 $G4VIS_BUILD_OPACS_DRIVER
setenv G4VIS_USE_OPACS	 $G4VIS_USE_OPACS
setenv G4VIS_BUILD_OPENGL_DRIVER	 $G4VIS_BUILD_OPENGL_DRIVER
setenv G4VIS_USE_OPENGL	 $G4VIS_USE_OPENGL
setenv G4VIS_BUILD_OPENGLX_DRIVER	 $G4VIS_BUILD_OPENGLX_DRIVER
setenv G4VIS_USE_OPENGLX	 $G4VIS_USE_OPENGLX
setenv G4VIS_BUILD_OIX_DRIVER	 $G4VIS_BUILD_OIX_DRIVER
setenv G4VIS_USE_OIX	 $G4VIS_USE_OIX
setenv G4VIS_BUILD_RAYTRACER_DRIVER	 $G4VIS_BUILD_RAYTRACER_DRIVER
setenv G4VIS_USE_RAYTRACER	 $G4VIS_USE_RAYTRACER
setenv G4VIS_BUILD_VRML_DRIVER	 $G4VIS_BUILD_VRML_DRIVER
setenv G4VIS_USE_VRML	 $G4VIS_USE_VRML
setenv G4VIS_BUILD_VRMLFILE_DRIVER	 $G4VIS_BUILD_VRMLFILE_DRIVER
setenv G4VIS_USE_VRMLFILE	 $G4VIS_USE_VRMLFILE
setenv G4VIS_BUILD_ASCIITREE_DRIVER	 $G4VIS_BUILD_ASCIITREE_DRIVER
setenv G4VIS_USE_ASCIITREE	 $G4VIS_USE_ASCIITREE

setenv G4LIB_BUILD_SHARED	 $G4LIB_BUILD_SHARED
setenv G4LIB_BUILD_STATIC	 $G4LIB_BUILD_STATIC
setenv G4LIB_USE_SHARED	 $G4LIB_USE_SHARED
setenv G4LIB_USE_STATIC	 $G4LIB_USE_STATIC

###############################################


#####################################################################
echo "On this machine the G4SYSTEM=\$G4SYSTEM"

echo "On this machine the G4INSTALL=\$G4INSTALL"
echo "On this machine the G4INCLUDE=\$G4INCLUDE"
echo "On this machine the G4BASE=\$G4BASE"
echo "On this machine the G4WORKDIR=\$G4WORKDIR"
echo "On this machine the G4TMP=\$G4TMP"
echo "On this machine the G4LIB=\$G4LIB"
echo "On this machine the G4BIN=\$G4BIN"
echo "On this machine the G4DATA=\$G4DATA"
echo "On this machine the G4LEVELGAMMADATA=\$G4LEVELGAMMADATA"

echo "On this machine the CLHEP_BASE_DIR=\$CLHEP_BASE_DIR"
echo "On this machine the CLHEP_INCLUDE_DIR=\$CLHEP_INCLUDE_DIR"
echo "On this machine the CLHEP_LIB_DIR=\$CLHEP_LIB_DIR"
echo "On this machine the CLHEP_LIB=\$CLHEP_LIB"

echo "On this machine the G4USE_OSPACE=\$G4USE_OSPACE"
echo "On this machine the G4OSPACE_BASE_DIR=\$G4OSPACE_BASE_DIR"

echo "On this machine the G4LIB_USE_GRANULAR=\$G4LIB_USE_GRANULAR"

echo "On this machine the G4DEBUG=\$G4DEBUG"

echo "On this machine the G4UI_BUILD_TERMINAL_SESSION=\$G4UI_BUILD_TERMINAL_SESSION"
echo "On this machine the G4UI_USE_TERMINAL=\$G4UI_USE_TERMINAL"
echo "On this machine the G4UI_BUILD_GAG_SESSION=\$G4UI_BUILD_GAG_SESSION"
echo "On this machine the G4UI_USE_GAG=\$G4UI_USE_GAG"
echo "On this machine the G4UI_BUILD_XAW_SESSION=\$G4UI_BUILD_XAW_SESSION"
echo "On this machine the G4UI_USE_XAW=\$G4UI_USE_XAW"
echo "On this machine the G4UI_BUILD_XM_SESSION=\$G4UI_BUILD_XM_SESSION"
echo "On this machine the G4UI_USE_XM=\$G4UI_USE_XM"
echo "On this machine the G4UI_BUILD_WIN32_SESSION=\$G4UI_BUILD_WIN32_SESSION"
echo "On this machine the G4UI_USE_WIN32=\$G4UI_USE_WIN32"
echo "On this machine the G4UI_USE_TCSH=\$G4UI_USE_TCSH"

echo "On this machine the G4ANALYSIS_BUILD=\$G4ANALYSIS_BUILD"
echo "On this machine the G4ANALYSIS_BUILD_LAB=\$G4ANALYSIS_BUILD_LAB"
echo "On this machine the G4ANALYSIS_BUILD_JAS=\$G4ANALYSIS_BUILD_JAS"
echo "On this machine the G4ANALYSIS_BUILD_LIZARD=\$G4ANALYSIS_BUILD_LIZARD"
echo "On this machine the G4ANALYSIS_USE=\$G4ANALYSIS_USE"
echo "On this machine the G4ANALYSIS_USE_LAB=\$G4ANALYSIS_USE_LAB"
echo "On this machine the G4ANALYSIS_USE_JAS=\$G4ANALYSIS_USE_JAS"
echo "On this machine the G4ANALYSIS_USE_LIZARD=\$G4ANALYSIS_USE_LIZARD"

echo "On this machine the G4USE_HEPODBMS=\$G4USE_HEPODBMS"

echo "On this machine the G4VIS_BUILD_DAWN_DRIVER=\$G4VIS_BUILD_DAWN_DRIVER"
echo "On this machine the G4VIS_USE_DAWN=\$G4VIS_USE_DAWN"
echo "On this machine the G4VIS_BUILD_DAWNFILE_DRIVER=\$G4VIS_BUILD_DAWNFILE_DRIVER"
echo "On this machine the G4VIS_USE_DAWNFILE=\$G4VIS_USE_DAWNFILE"
echo "On this machine the G4VIS_BUILD_OPACS_DRIVER=\$G4VIS_BUILD_OPACS_DRIVER"
echo "On this machine the G4VIS_USE_OPACS=\$G4VIS_USE_OPACS"
echo "On this machine the G4VIS_BUILD_OPENGL_DRIVER=\$G4VIS_BUILD_OPENGL_DRIVER"
echo "On this machine the G4VIS_USE_OPENGL=\$G4VIS_USE_OPENGL"
echo "On this machine the G4VIS_BUILD_OPENGLX_DRIVER=\$G4VIS_BUILD_OPENGLX_DRIVER"
echo "On this machine the G4VIS_USE_OPENGLX=\$G4VIS_USE_OPENGLX"
echo "On this machine the G4VIS_BUILD_OIX_DRIVER=\$G4VIS_BUILD_OIX_DRIVER"
echo "On this machine the G4VIS_USE_OIX=\$G4VIS_USE_OIX"
echo "On this machine the G4VIS_BUILD_RAYTRACER_DRIVER=\$G4VIS_BUILD_RAYTRACER_DRIVER"
echo "On this machine the G4VIS_USE_RAYTRACER=\$G4VIS_USE_RAYTRACER"
echo "On this machine the G4VIS_BUILD_VRML_DRIVER=\$G4VIS_BUILD_VRML_DRIVER"
echo "On this machine the G4VIS_USE_VRML=\$G4VIS_USE_VRML"
echo "On this machine the G4VIS_BUILD_VRMLFILE_DRIVER=\$G4VIS_BUILD_VRMLFILE_DRIVER"
echo "On this machine the G4VIS_USE_VRMLFILE=\$G4VIS_USE_VRMLFILE"
echo "On this machine the G4VIS_BUILD_ASCIITREE_DRIVER=\$G4VIS_BUILD_ASCIITREE_DRIVER"
echo "On this machine the G4VIS_USE_ASCIITREE=\$G4VIS_USE_ASCIITREE"

echo "On this machine the G4LIB_BUILD_SHARED=\$G4LIB_BUILD_SHARED"
echo "On this machine the G4LIB_BUILD_STATIC=\$G4LIB_BUILD_STATIC"
echo "On this machine the G4LIB_USE_SHARED=\$G4LIB_USE_SHARED"
echo "On this machine the G4LIB_USE_STATIC=\$G4LIB_USE_STATIC"

###############################################################
!GROK!THIS!

: In the following dollars and backticks do not need the extra backslash.
$spitshell >>env.csh <<'!NO!SUBS!'
!NO!SUBS!
chmod 755 env.csh
$eunicefix env.csh
