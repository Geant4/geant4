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
G4SYSTEM="$g4system"

#
# g4dirs.U
#
G4INSTALL="$g4install"
G4INCLUDE="$g4include"
G4BASE="$g4base"
G4WORKDIR="$g4workdir"
G4TMP="$g4tmp"
G4LIB="$g4lib"
G4BIN="$g4bin"
G4DATA="$g4data"
G4LEVELGAMMADATA="$g4levelgammadata"

#
# g4clhep.U
#
CLHEP_BASE_DIR="$g4clhep_base_dir"
CLHEP_INCLUDE_DIR="$g4clhep_include_dir"
CLHEP_LIB_DIR="$g4clhep_lib_dir"
CLHEP_LIB="$g4clhep_lib"

#
# g4ospace
#
G4USE_OSPACE="$g4use_ospace"
G4OSPACE_BASE_DIR="$g4ospace_base_dir"

#
# g4debug
#
G4DEBUG="$g4debug"

#
# g4analysis
#
if [ X$g4analysis_build = Xy ] ; then
G4ANALYSIS_BUILD=1
else
G4ANALYSIS_BUILD=""
fi 

if [ X$g4analysis_build_lab = Xy ] ; then
G4ANALYSIS_BUILD_LAB=1
else
G4ANALYSIS_BUILD_LAB=""
fi 

if [ X$g4analysis_build_jas = Xy ] ; then
G4ANALYSIS_BUILD_JAS=1
else
G4ANALYSIS_BUILD_JAS=""
fi 

if [ X$g4analysis_build_lizard = Xy ] ; then
G4ANALYSIS_BUILD_LIZARD=1
else
G4ANALYSIS_BUILD_LIZARD=""
fi 

if [ X$g4analysis_use = Xy ] ; then
G4ANALYSIS_USE=1
else
G4ANALYSIS_USE=""
fi 

if [ X$g4analysis_use_lab = Xy ] ; then
G4ANALYSIS_USE_LAB=1
else
G4ANALYSIS_USE_LAB=""
fi 

if [ X$g4analysis_use_jas = Xy ] ; then
G4ANALYSIS_USE_JAS=1
else
G4ANALYSIS_USE_JAS=""
fi 

if [ X$g4analysis_use_lizard = Xy ] ; then
G4ANALYSIS_USE_LIZARD=1
else
G4ANALYSIS_USE_LIZARD=""
fi 

#
# g4ui
#
if [ X$g4ui_build_terminal_session = Xy ] ; then
G4UI_BUILD_TERMINAL_SESSION=1
else
G4UI_BUILD_TERMINAL_SESSION=""
fi 

if [ X$g4ui_use_terminal = Xy ] ; then
G4UI_USE_TERMINAL=1
else
G4UI_USE_TERMINAL=""
fi 

if [ X$g4ui_build_gag_session = Xy ] ; then
G4UI_BUILD_GAG_SESSION=1
else
G4UI_BUILD_GAG_SESSION=""
fi 

if [ X$g4ui_use_gag = Xy ] ; then
G4UI_USE_GAG=1
else
G4UI_USE_GAG=""
fi 

if [ X$g4ui_build_xaw_session = Xy ] ; then
G4UI_BUILD_XAW_SESSION=1
else
G4UI_BUILD_XAW_SESSION=""
fi 

if [ X$g4ui_use_xaw = Xy ] ; then
G4UI_USE_XAW=1
else
G4UI_USE_XAW=""
fi 

if [ X$g4ui_build_xm_session = Xy ] ; then
G4UI_BUILD_XM_SESSION=1
else
G4UI_BUILD_XM_SESSION=""
fi 

if [ X$g4ui_use_xm = Xy ] ; then
G4UI_USE_XM=1
else
G4UI_USE_XM=""
fi 

if [ X$g4ui_build_win32_session = Xy ] ; then
G4UI_BUILD_WIN32_SESSION=1
else
G4UI_BUILD_WIN32_session=""
fi 

if [ X$g4ui_use_win32 = Xy ] ; then
G4UI_USE_WIN32=1
else
G4UI_USE_WIN32=""
fi 

if [ X$g4ui_use_tcsh = Xy ] ; then
G4UI_USE_TCSH=1
else
G4UI_USE_TCSH=""
fi 
#
# g4shared
#
G4LIB_BUILD_SHARED="$g4lib_build_shared"
G4LIB_BUILD_STATIC="$g4lib_build_static"
G4LIB_USE_SHARED="$g4lib_use_shared"
G4LIB_USE_STATIC="$g4lib_use_static"

if [ X$g4lib_build_shared = Xy ] ; then
G4LIB_BUILD_SHARED=1
else
G4LIB_BUILD_SHARED=""
fi 

if [ X$g4lib_build_static = Xy ] ; then
G4LIB_BUILD_STATIC=1
else
G4LIB_BUILD_STATIC=""
fi 

if [ X$g4lib_use_shared = Xy ] ; then
G4LIB_USE_SHARED=1
else
G4LIB_USE_SHARED=""
fi 

if [ X$g4lib_use_static = Xy ] ; then
G4LIB_USE_STATIC=1
else
G4LIB_USE_STATIC=""
fi 

#
# g4granular
#
if [ X$g4lib_use_granular = Xy ] ; then
G4LIB_USE_GRANULAR=1
else
G4LIB_USE_GRANULAR=""
fi 

#####################################################################
export G4SYSTEM

export G4INSTALL
export G4INCLUDE
export G4BASE
export G4WORKDIR
export G4TMP
export G4LIB
export G4BIN
export G4DATA
export G4LEVELGAMMADATA

export CLHEP_BASE_DIR
export CLHEP_INCLUDE_DIR
export CLHEP_LIB_DIR
export CLHEP_LIB

export G4USE_OSPACE
export G4OSPACE_BASE_DIR

export G4LIB_USE_GRANULAR

export G4DEBUG

export G4ANALYSIS_BUILD
export G4ANALYSIS_BUILD_LAB
export G4ANALYSIS_BUILD_JAS
export G4ANALYSIS_BUILD_LIZARD
export G4ANALYSIS_USE
export G4ANALYSIS_USE_LAB
export G4ANALYSIS_USE_JAS
export G4ANALYSIS_USE_LIZARD

export G4UI_BUILD_TERMINAL_SESSION
export G4UI_USE_TERMINAL
export G4UI_BUILD_GAG_SESSION
export G4UI_USE_GAG
export G4UI_BUILD_XAW_SESSION
export G4UI_USE_XAW
export G4UI_BUILD_XM_SESSION
export G4UI_USE_XM
export G4UI_BUILD_WIN32_SESSION
export G4UI_USE_WIN32
export G4UI_USE_TCSH

export G4LIB_BUILD_SHARED
export G4LIB_BUILD_STATIC
export G4LIB_USE_SHARED
export G4LIB_USE_STATIC

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

echo "On this machine the CLHEP_BASE_DIR=\$G4CLHEP_BASE_DIR"
echo "On this machine the CLHEP_INCLUDE_DIR=\$G4CLHEP_INCLUDE_DIR"
echo "On this machine the CLHEP_LIB_DIR=\$G4CLHEP_LIB_DIR"
echo "On this machine the CLHEP_LIB=\$G4CLHEP_LIB"

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

echo "On this machine the G4LIB_BUILD_SHARED=\$G4LIB_BUILD_SHARED"
echo "On this machine the G4LIB_BUILD_STATIC=\$G4LIB_BUILD_STATIC"
echo "On this machine the G4LIB_USE_SHARED=\$G4LIB_USE_SHARED"
echo "On this machine the G4LIB_USE_STATIC=\$G4LIB_USE_STATIC"

###############################################################
echo ""
echo "OK, going to $g4install/source and start 'gmake'..."
echo ""
cd $g4install/source
echo `pwd`
if [ X$g4global = Xy ] ; then
gmake global
if [ X$g4granular = Xy ] ; then
gmake
fi
else
gmake
fi
if [ X$g4includes_flag = Xy ] ; then
gmake includes
fi
!GROK!THIS!

: In the following dollars and backticks do not need the extra backslash.
$spitshell >>install.sh <<'!NO!SUBS!'
!NO!SUBS!
chmod 755 install.sh
$eunicefix install.sh
