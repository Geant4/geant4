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
echo "Extracting install.csh (with variable substitutions)"
: This section of the file will have variable substitutions done on it.
: Move anything that needs config subs from !NO!SUBS! section to !GROK!THIS!.
: Protect any dollar signs and backticks that you do not want interpreted
: by putting a backslash in front.  You may delete these comments.
$spitshell >install.csh <<!GROK!THIS!
#
# g4system.U
#
set g4system="$g4system"

#
# g4dirs.U
#
set g4install="$g4install"
set g4base="$g4base"
set g4workdir="$g4workdir"
set g4tmp="$g4tmp"
set g4lib="$g4lib"
set g4bin="$g4bin"
set g4data="$g4data"
set g4levelgammadata="$g4levelgammadata"

#
# g4clhep.U
#
set g4clhep_base_dir="$g4clhep_base_dir"
set g4clhep_include_dir="$g4clhep_include_dir"
set g4clhep_lib_dir="$g4clhep_lib_dir"
set g4clhep_lib="$g4clhep_lib"

#
# g4ospace.U
#
set g4use_ospace="$g4use_ospace"
set g4ospace_base_dir="$g4ospace_base_dir"

#################################################
setenv G4SYSTEM "$g4system"

setenv G4INSTALL "$g4install"
setenv G4BASE "$g4base"
setenv G4WORKDIR "$g4workdir"
setenv G4TMP "$g4tmp"
setenv G4LIB "$g4lib"
setenv G4BIN "$g4bin"
setenv G4DATA "$g4data"
setenv G4LEVELGAMMADATA "$g4levelgammadata"

setenv G4CLHEP_BASE_DIR "$g4clhep_base_dir"
setenv G4CLHEP_INCLUDE_DIR "$g4clhep_include_dir"
setenv G4CLHEP_LIB_DIR "$g4clhep_lib_dir"
setenv G4CLHEP_LIB "$g4clhep_lib"

setenv G4USE_OSPACE "$g4use_ospace"
setenv G4OSPACE_BASE_DIR "$g4ospace_base_dir"

####################################################
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
echo ""
echo "OK, going to $g4install/source and start 'make'..."
echo ""
cd $g4install/source
echo `pwd`
if [ X$g4global = Xglobal ] ; then
make global
else
make
fi
!GROK!THIS!

: In the following dollars and backticks do not need the extra backslash.
$spitshell >>install.csh <<'!NO!SUBS!'
!NO!SUBS!
chmod 755 install.csh
$eunicefix install.csh
