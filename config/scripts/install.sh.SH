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
echo "Extracting install.sh (with variable substitutions)"
: This section of the file will have variable substitutions done on it.
: Move anything that needs config subs from !NO!SUBS! section to !GROK!THIS!.
: Protect any dollar signs and backticks that you do not want interpreted
: by putting a backslash in front.  You may delete these comments.
$spitshell >install.sh <<!GROK!THIS!
#
# g4system.U
#
g4system="$g4system"

#
# g4dirs.U
#
g4install="$g4install"
g4include="$g4include"
g4base="$g4base"
g4workdir="$g4workdir"
g4tmp="$g4tmp"
g4lib="$g4lib"
g4bin="$g4bin"
g4data="$g4data"
g4levelgammadata="$g4levelgammadata"

#
# g4clhep.U
#
g4clhep_base_dir="$g4clhep_base_dir"
g4clhep_include_dir="$g4clhep_include_dir"
g4clhep_lib_dir="$g4clhep_lib_dir"
g4clhep_lib="$g4clhep_lib"

#
# g4ospace
#
g4use_ospace="$g4use_ospace"
g4ospace_base_dir="$g4ospace_base_dir"

#
# g4granular
#
g4global="$g4global"
g4lib_use_granular="$g4lib_use_granular"

#
# g4debug
#
g4debug="$g4debug"

#
# g4shared
#
if [ X$g4lib_build_shared = Xy ] ; then
g4lib_build_shared=1
else
g4lib_build_shared=""
fi 

if [ X$g4lib_build_static = Xy ] ; then
g4lib_build_static=1
else
g4lib_build_static=""
fi 

if [ X$g4lib_use_shared = Xy ] ; then
g4lib_use_shared=1
else
g4lib_use_shared=""
fi 

if [ X$g4lib_use_static = Xy ] ; then
g4lib_use_static=1
else
g4lib_use_static=""
fi 

!GROK!THIS!

: In the following dollars and backticks do not need the extra backslash.
$spitshell >>install.sh <<'!NO!SUBS!'
######################################
G4SYSTEM="$g4system"

G4INSTALL="$g4install"
G4INCLUDE="$g4include"
G4BASE="$g4base"
G4WORKDIR="$g4workdir"
G4TMP="$g4tmp"
G4LIB="$g4lib"
G4BIN="$g4bin"
G4DATA="$g4data"
G4LEVELGAMMADATA="$g4levelgammadata"

G4CLHEP_BASE_DIR="$g4clhep_base_dir"
G4CLHEP_INCLUDE_DIR="$g4clhep_include_dir"
G4CLHEP_LIB_DIR="$g4clhep_lib_dir"
G4CLHEP_LIB="$g4clhep_lib"

G4USE_OSPACE="$g4use_ospace"
G4OSPACE_BASE_DIR="$g4ospace_base_dir"

G4LIB_USE_GRANULAR="$g4lib_use_granular"

G4DEBUG="$g4debug"

G4LIB_BUILD_SHARED="$g4lib_build_shared"
G4LIB_BUILD_STATIC="$g4lib_build_static"
G4LIB_USE_SHARED="$g4lib_use_shared"
G4LIB_USE_STATIC="$g4lib_use_static"

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

export G4CLHEP_BASE_DIR
export G4CLHEP_INCLUDE_DIR
export G4CLHEP_LIB_DIR
export G4CLHEP_LIB

export G4USE_OSPACE
export G4OSPACE_BASE_DIR

export G4LIB_USE_GRANULAR

export G4DEBUG

export G4LIB_BUILD_SHARED
export G4LIB_BUILD_STATIC
export G4LIB_USE_SHARED
export G4LIB_USE_STATIC

#####################################################################
echo "On this machine the G4SYSTEM=$g4system"

echo "On this machine the G4INSTALL=$g4install"
echo "On this machine the G4INCLUDE=$g4include"
echo "On this machine the G4BASE=$g4base"
echo "On this machine the G4WORKDIR=$g4workdir"
echo "On this machine the G4TMP=$g4tmp"
echo "On this machine the G4LIB=$g4lib"
echo "On this machine the G4BIN=$g4bin"
echo "On this machine the G4DATA=$g4data"
echo "On this machine the G4LEVELGAMMADATA=$g4levelgammadata"

echo "On this machine the G4CLHEP_BASE_DIR=$g4clhep_base_dir"
echo "On this machine the G4CLHEP_INCLUDE_DIR=$g4clhep_include_dir"
echo "On this machine the G4CLHEP_LIB_DIR=$g4clhep_lib_dir"
echo "On this machine the G4CLHEP_LIB=$g4clhep_lib"

echo "On this machine the G4USE_OSPACE=$g4use_ospace"
echo "On this machine the G4OSPACE_BASE_DIR=$g4ospace_base_dir"

echo "On this machine the G4LIB_USE_GRANULAR=$g4lib_use_granular"

echo "On this machine the G4DEBUG=$g4debug"

echo "On this machine the G4LIB_BUILD_SHARED=$g4lib_build_shared"
echo "On this machine the G4LIB_BUILD_STATIC=$g4lib_build_static"
echo "On this machine the G4LIB_USE_SHARED=$g4lib_use_shared"
echo "On this machine the G4LIB_USE_STATIC=$g4lib_use_static"

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
if [ X$g4include = X1 ] ; then
gmake includes
fi
!NO!SUBS!
chmod 755 install.sh
$eunicefix install.sh
