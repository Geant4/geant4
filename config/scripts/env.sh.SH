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
#echo "Extracting env.sh (with variable substitutions)"
: This section of the file will have variable substitutions done on it.
: Move anything that needs config subs from !NO!SUBS! section to !GROK!THIS!.
: Protect any dollar signs and backticks that you do not want interpreted
: by putting a backslash in front.  You may delete these comments.
$spitshell >env.sh <<!GROK!THIS!
#
# g4system.U
#
g4system="$g4system"

#
# g4dirs.U
#
g4install="$g4install"
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

######################################
export G4SYSTEM="$g4system"

export G4INSTALL="$g4install"
export G4BASE="$g4base"
export G4WORKDIR="$g4workdir"
export G4TMP="$g4tmp"
export G4LIB="$g4lib"
export G4BIN="$g4bin"
export G4DATA="$g4data"
export G4LEVELGAMMADATA="$g4levelgammadata"

export G4CLHEP_BASE_DIR="$g4clhep_base_dir"
export G4CLHEP_INCLUDE_DIR="$g4clhep_include_dir"
export G4CLHEP_LIB_DIR="$g4clhep_lib_dir"
export G4CLHEP_LIB="$g4clhep_lib"

export G4USE_OSPACE="$g4use_ospace"
export G4OSPACE_BASE_DIR="$g4ospace_base_dir"

#####################################################################
echo "On this machine the G4SYSTEM=$g4system"

echo "On this machine the G4INSTALL=$g4install"
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
###############################################################
!GROK!THIS!

: In the following dollars and backticks do not need the extra backslash.
$spitshell >>env.sh <<'!NO!SUBS!'
!NO!SUBS!
chmod 755 env.sh
$eunicefix env.sh
