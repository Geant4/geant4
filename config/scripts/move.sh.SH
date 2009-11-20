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
#echo "Extracting move.sh (with variable substitutions)"
: This section of the file will have variable substitutions done on it.
: Move anything that needs config subs from !NO!SUBS! section to !GROK!THIS!.
: Protect any dollar signs and backticks that you do not want interpreted
: by putting a backslash in front.  You may delete these comments.
$spitshell >move.sh <<!GROK!THIS!
######################################
#
# g4system.U
#
#+
G4SYSTEM="$g4system"
export G4SYSTEM
echo "On this machine the G4SYSTEM=\$G4SYSTEM"

#
# g4dirs.U
#
#+
G4INSTALL="$g4install"
export G4INSTALL
echo "On this machine the G4INSTALL=\$G4INSTALL"

#+
if [ X$g4include != X -a X$g4include != Xn ] ; then
G4INCLUDE="$g4include"
export G4INCLUDE
echo "On this machine the G4INCLUDE=\$G4INCLUDE"
fi

#+
if [ X$g4tmp != X ] ; then
G4TMP="$g4tmp"
export G4TMP
echo "On this machine the G4TMP=\$G4TMP"
fi

#+
if [ X$g4lib != X ] ; then 
G4LIB="$g4lib"
export G4LIB
echo "On this machine the G4LIB=\$G4LIB"
fi

#+
if [ X$g4lib != X ] ; then 
G4BIN="$g4lib/../bin"
export G4BIN
fi

#+
if [ X$g4levelgammadata != X ] ; then 
G4LEVELGAMMADATA="$g4levelgammadata"
export G4LEVELGAMMADATA
echo "On this machine the G4LEVELGAMMADATA=\$G4LEVELGAMMADATA"
fi

#+
if [ X$g4radioactivedata != X ] ; then 
G4RADIOACTIVEDATA="$g4radioactivedata"
export G4RADIOACTIVEDATA
echo "On this machine the G4RADIOACTIVEDATA=\$G4RADIOACTIVEDATA"
fi

#+
if [ X$g4ledata != X ] ; then 
G4LEDATA="$g4ledata"
export G4LEDATA
echo "On this machine the G4LEDATA=\$G4LEDATA"
fi

#+
if [ X$g4neutronhpcrosssections != X ] ; then 
G4NEUTRONHPDATA="$g4neutronhpcrosssections"
export G4NEUTRONHPDATA
echo "On this machine the G4NEUTRONHPDATA=\$G4NEUTRONHPDATA"
fi

#+
if [ X$g4abladata != X ] ; then 
G4ABLADATA="$g4abladata"
export G4ABLADATA
echo "On this machine the G4ABLADATA=\$G4ABLADATA"
fi



#+
#if [ X$g4elasticdata != X ] ; then 
#G4ELASTICDATA="$g4elasticdata"
#export G4ELASTICDATA
#echo "On this machine the G4ELASTICDATA=\$G4ELASTICDATA"
#fi

#
# g4clhep.U
#
if [ X$g4clhep_base_dir != X/usr -o X$g4clhep_base_dir != X/usr/local ] ; then 
CLHEP_BASE_DIR="$g4clhep_base_dir"
export CLHEP_BASE_DIR
echo "On this machine the CLHEP_BASE_DIR=\$CLHEP_BASE_DIR"
fi

#+
if [ X$g4clhep_include_dir != X ] ; then 
CLHEP_INCLUDE_DIR="$g4clhep_include_dir"
export CLHEP_INCLUDE_DIR
echo "On this machine the CLHEP_INCLUDE_DIR=\$CLHEP_INCLUDE_DIR"
fi

#+
if [ X$g4clhep_lib_dir != X ] ; then 
CLHEP_LIB_DIR="$g4clhep_lib_dir"
export CLHEP_LIB_DIR
echo "On this machine the CLHEP_LIB_DIR=\$CLHEP_LIB_DIR"
fi

#+
if [ X$g4clhep_lib != X ] ; then 
CLHEP_LIB="$g4clhep_lib"
export CLHEP_LIB
echo "On this machine the CLHEP_LIB=\$CLHEP_LIB"
fi

#+
#
# g4debug
#
if [ X$g4debug = Xy ] ; then 
G4DEBUG=1
export G4DEBUG
echo "On this machine the G4DEBUG=\$G4DEBUG"
fi


#
# g4ui
#
#+
if [ X$g4ui_build_xaw_session = Xy ] ; then
G4UI_BUILD_XAW_SESSION=1
export G4UI_BUILD_XAW_SESSION
echo "On this machine the G4UI_BUILD_XAW_SESSION=\$G4UI_BUILD_XAW_SESSION"
fi 

#+
if [ X$g4ui_use_xaw = Xy ] ; then
G4UI_USE_XAW=1
export G4UI_USE_XAW
echo "On this machine the G4UI_USE_XAW=\$G4UI_USE_XAW"
fi 

#+
if [ X$g4ui_build_xm_session = Xy ] ; then
G4UI_BUILD_XM_SESSION=1
export G4UI_BUILD_XM_SESSION
echo "On this machine the G4UI_BUILD_XM_SESSION=\$G4UI_BUILD_XM_SESSION"
fi 

#+
if [ X$g4ui_use_xm = Xy ] ; then
G4UI_USE_XM=1
export G4UI_USE_XM
echo "On this machine the G4UI_USE_XM=\$G4UI_USE_XM"
fi 

#+
if [ X$g4ui_build_win32_session = Xy ] ; then
G4UI_BUILD_WIN32_SESSION=1
export G4UI_BUILD_WIN32_SESSION
echo "On this machine the G4UI_BUILD_WIN32_SESSION=\$G4UI_BUILD_WIN32_SESSION"
fi 

#+
if [ X$g4ui_use_win32 = Xy ] ; then
G4UI_USE_WIN32=1
export G4UI_USE_WIN32
echo "On this machine the G4UI_USE_WIN32=\$G4UI_USE_WIN32"
fi 


#+
if [ X$g4ui_build_qt_session = Xy ] ; then
G4UI_BUILD_QT_SESSION=1
export G4UI_BUILD_QT_SESSION
echo "On this machine the G4UI_BUILD_QT_SESSION=\$G4UI_BUILD_QT_SESSION"
fi 

#+
if [ X$g4ui_use_qt = Xy ] ; then
G4UI_USE_QT=1
export G4UI_USE_QT
echo "On this machine the G4UI_USE_QT=\$G4UI_USE_QT"
fi 



#
# g4vis
#
#+
if [ X$g4vis_build_dawn_driver = Xy ] ; then
G4VIS_BUILD_DAWN_DRIVER=1
export G4VIS_BUILD_DAWN_DRIVER
echo "On this machine the G4VIS_BUILD_DAWN_DRIVER=\$G4VIS_BUILD_DAWN_DRIVER"
fi 

#+
if [ X$g4vis_build_openglx_driver = Xy ] ; then
G4VIS_BUILD_OPENGLX_DRIVER=1
export G4VIS_BUILD_OPENGLX_DRIVER
echo "On this machine the G4VIS_BUILD_OPENGLX_DRIVER=\$G4VIS_BUILD_OPENGLX_DRIVER"
fi 

#+
if [ X$g4vis_build_openglxm_driver = Xy ] ; then
G4VIS_BUILD_OPENGLXM_DRIVER=1
export G4VIS_BUILD_OPENGLXM_DRIVER
echo "On this machine the G4VIS_BUILD_OPENGLXM_DRIVER=\$G4VIS_BUILD_OPENGLXM_DRIVER"
fi 

#+
if [ X$g4vis_build_openglwin32_driver = Xy ] ; then
G4VIS_BUILD_OPENGLWIN32_DRIVER=1
export G4VIS_BUILD_OPENGLWIN32_DRIVER
echo "On this machine the G4VIS_BUILD_OPENGLWIN32_DRIVER=\$G4VIS_BUILD_OPENGLWIN32_DRIVER"
fi 

#+
if [ X$g4vis_build_oix_driver = Xy ] ; then
G4VIS_BUILD_OIX_DRIVER=1
export G4VIS_BUILD_OIX_DRIVER
echo "On this machine the G4VIS_BUILD_OIX_DRIVER=\$G4VIS_BUILD_OIX_DRIVER"
fi 

#+
if [ X$g4vis_build_oiwin32_driver = Xy ] ; then
G4VIS_BUILD_OIWIN32_DRIVER=1
export G4VIS_BUILD_OIWIN32_DRIVER
echo "On this machine the G4VIS_BUILD_OIWIN32_DRIVER=\$G4VIS_BUILD_OIWIN32_DRIVER"
fi 

#+
if [ X$g4vis_build_vrml_driver = Xy ] ; then
G4VIS_BUILD_VRML_DRIVER=1
export G4VIS_BUILD_VRML_DRIVER
echo "On this machine the G4VIS_BUILD_VRML_DRIVER=\$G4VIS_BUILD_VRML_DRIVER"
fi 

#+
if [ X$g4vis_build_openglqt_driver = Xy ] ; then
G4VIS_BUILD_OPENGLQT_DRIVER=1
export G4VIS_BUILD_OPENGLQT_DRIVER
echo "On this machine the G4VIS_BUILD_OPENGLQT_DRIVER=\$G4VIS_BUILD_OPENGLQT_DRIVER"
fi 



#+
if [ X$g4vis_use_dawn = Xy ] ; then
G4VIS_USE_DAWN=1
export G4VIS_USE_DAWN
echo "On this machine the G4VIS_USE_DAWN=\$G4VIS_USE_DAWN"
fi 

#+
if [ X$g4vis_use_openglx = Xy ] ; then
G4VIS_USE_OPENGLX=1
export G4VIS_USE_OPENGLX
echo "On this machine the G4VIS_USE_OPENGLX=\$G4VIS_USE_OPENGLX"
fi 

#+
if [ X$g4vis_use_openglxm = Xy ] ; then
G4VIS_USE_OPENGLXM=1
export G4VIS_USE_OPENGLXM
echo "On this machine the G4VIS_USE_OPENGLXM=\$G4VIS_USE_OPENGLXM"
fi 

#+
if [ X$g4vis_use_openglwin32 = Xy ] ; then
G4VIS_USE_OPENGLWIN32=1
export G4VIS_USE_OPENGLWIN32
echo "On this machine the G4VIS_USE_OPENGLWIN32=\$G4VIS_USE_OPENGLWIN32"
fi 

#+
if [ X$g4vis_use_oix = Xy ] ; then
G4VIS_USE_OIX=1
export G4VIS_USE_OIX
echo "On this machine the G4VIS_USE_OIX=\$G4VIS_USE_OIX"
fi 

#+
if [ X$g4vis_use_oiwin32 = Xy ] ; then
G4VIS_USE_OIWIN32=1
export G4VIS_USE_OIWIN32
echo "On this machine the G4VIS_USE_OIWIN32=\$G4VIS_USE_OIWIN32"
fi 

#+
if [ X$g4vis_use_vrml = Xy ] ; then
G4VIS_USE_VRML=1
export G4VIS_USE_VRML
echo "On this machine the G4VIS_USE_VRML=\$G4VIS_USE_VRML"
fi 

#+
if [ X$g4vis_use_openglqt = Xy ] ; then
G4VIS_USE_OPENGLQT=1
export G4VIS_USE_OPENGLQT
echo "On this machine the G4VIS_USE_OPENGLQT=\$G4VIS_USE_OPENGLQT"
fi 



#+
if [ X$g4vis_oglhome != X ] ; then
OGLHOME="$g4vis_oglhome"
export OGLHOME
echo "On this machine the OGLHOME=\$OGLHOME"
fi 

#+
if [ X$g4_qthome != X ] ; then
QTHOME=$g4_qthome
export QTHOME
echo "On this machine the QTHOME=\$QTHOME"
fi 



#
# Use GDML module
#
#+
if [ X$g4lib_build_gdml = Xy ] ; then
G4LIB_BUILD_GDML=1
export G4LIB_BUILD_GDML
echo "On this machine the G4LIB_BUILD_GDML=\$G4LIB_BUILD_GDML"
fi 

if [ X$g4lib_build_gdml = Xy ] ; then
XERCESCROOT=$g4gdml_xercesc_root
export XERCESCROOT
echo "On this machine the XERCESCROOT=\$XERCESCROOT"
fi 



#
# Use G3TOG4 module
#
#+
if [ X$g4wlib_build_g3tog4 = Xy ] ; then
G4LIB_BUILD_G3TOG4=1
export G4LIB_BUILD_G3TOG4
echo "On this machine the G4LIB_BUILD_G3TOG4=\$G4LIB_BUILD_G3TOG4"
fi 

if [ X$g4w_use_g3tog4 = Xy ] ; then
G4LIB_USE_G3TOG4=1
export G4LIB_USE_G3TOG4
echo "On this machine the G4LIB_USE_G3TOG4=\$G4LIB_USE_G3TOG4"
fi 

#
# Use ZLIB module
#
#+
if [ X$g4wlib_build_zlib = Xy ] ; then
G4LIB_BUILD_ZLIB=1
export G4LIB_BUILD_ZLIB
echo "On this machine the G4LIB_BUILD_ZLIB=\$G4LIB_BUILD_ZLIB"
fi 

if [ X$g4w_use_zlib = Xy ] ; then
G4LIB_USE_ZLIB=1
export G4LIB_USE_ZLIB
echo "On this machine the G4LIB_USE_ZLIB=\$G4LIB_USE_ZLIB"
fi 


#+
#
# g4shared
#
if [ X$g4lib_build_shared = Xy ] ; then
G4LIB_BUILD_SHARED=1
export G4LIB_BUILD_SHARED
echo "On this machine the G4LIB_BUILD_SHARED=\$G4LIB_BUILD_SHARED"
fi 

if [ X$g4lib_build_static = Xy ] ; then
G4LIB_BUILD_STATIC=1
export G4LIB_BUILD_STATIC
echo "On this machine the G4LIB_BUILD_STATIC=\$G4LIB_BUILD_STATIC"
fi 

#+
#
# g4granular
#
if [ X$g4lib_use_granular = Xy ] ; then
G4LIB_USE_GRANULAR=1
export G4LIB_USE_GRANULAR
echo "On this machine the G4LIB_USE_GRANULAR=\$G4LIB_USE_GRANULAR"
fi 

#####################################################################

#####################################################################

echo ""
echo "Starting installation..."
echo ""

cd $g4install
#echo `pwd`


# Remove 'verbosity' in 'tar' command
# If not INCLUDE then install full tree in 'src' directory!
#
if [ X$g4includes_flag = Xn ] ; then
if [ X$g4install != X$g4final_install ] ; then
echo ""
echo "Installing sources and headers..."
echo ""
mkdir -p $g4final_install/src/geant4/
tar cf - Configure config/ source/ environments/ examples/ include/ |(cd $g4final_install/src/geant4/; tar xf -)
fi
fi

########################################
# There is no env!
if [ X$g4includes_flag = Xy ] ; then

if [ X$g4install != X$g4final_install ] ; then
echo ""
echo "Installing sources..."
echo ""

mkdir -p $g4final_install/src/geant4/
tar cf - Configure config/ source/ environments/ examples/ |(cd $g4final_install/src/geant4/; tar xf -)
fi

echo ""
echo "Installing headers..."
echo ""
cd $g4install/source
$g4make includes dependencies=""
fi

# Move libraries
if [ X$g4install != X$g4final_install ] ; then
echo ""
echo "Installing libraries..."
echo ""
mkdir -p $g4final_install/lib/geant4
cd $g4install/lib
tar cf - $g4system |(cd $g4final_install/lib/geant4; tar xf -)
fi

# For the reason of strange behaivior for DOTted directory...
mkdir -p $g4final_install/src/geant4/.config
cd $g4install/.config
tar cf - . |(cd $g4final_install/src/geant4/.config; tar xf -)


!GROK!THIS!

: In the following dollars and backticks do not need the extra backslash.
$spitshell >>move.sh <<'!NO!SUBS!'
!NO!SUBS!
chmod 755 move.sh
$eunicefix move.sh

