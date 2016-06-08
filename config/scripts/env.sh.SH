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
if [ X$g4include != X ] ; then 
G4INCLUDE="$g4include"
export G4INCLUDE
echo "On this machine the G4INCLUDE=\$G4INCLUDE"
fi

#+
if [ X$g4lib != X ] ; then 
G4LIB="$g4lib"
export G4LIB
echo "On this machine the G4LIB=\$G4LIB"
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
NeutronHPCrossSections="$g4neutronhpcrosssections"
export NeutronHPCrossSections
echo "On this machine the NeutronHPCrossSections=\$NeutronHPCrossSections"
fi

#
# g4clhep.U
#
if [ X$g4clhep_base_dir != X ] ; then 
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
if [ X$g4debug != X ] ; then 
G4DEBUG="$g4debug"
export G4DEBUG
echo "On this machine the G4DEBUG=\$G4DEBUG"
fi

#
# g4analysis
#
#+
if [ X$g4wanalysis_use = Xy ] ; then
G4ANALYSIS_USE=1
export G4ANALYSIS_USE
echo "On this machine the G4ANALYSIS_USE=\$G4ANALYSIS_USE"
fi 

#+
#
# g4persistency
#
if [ X$g4use_hepodbms = Xy ] ; then
G4USE_HEPODBMS=1
export G4USE_HEPODBMS
.  \$G4INSTALL/examples/extended/persistency/PersistentEx01/g4odbms_setup.sh
echo "On this machine the G4USE_HEPODBMS=\$G4USE_HEPODBMS"
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
if [ X$g4ui_build_wo_session = Xy ] ; then
G4UI_BUILD_WO_SESSION=1
export G4UI_BUILD_WO_SESSION
echo "On this machine the G4UI_BUILD_WO_SESSION=\$G4UI_BUILD_WO_SESSION"
fi 

#+
if [ X$g4ui_use_wo = Xy ] ; then
G4UI_USE_WO=1
export G4UI_USE_WO
echo "On this machine the G4UI_USE_WO=\$G4UI_USE_WO"
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
if [ X$g4vis_build_opacs_driver = Xy ] ; then
G4VIS_BUILD_OPACS_DRIVER=1
export G4VIS_BUILD_OPACS_DRIVER
echo "On this machine the G4VIS_BUILD_OPACS_DRIVER=\$G4VIS_BUILD_OPACS_DRIVER"
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
if [ X$g4vis_use_dawn = Xy ] ; then
G4VIS_USE_DAWN=1
export G4VIS_USE_DAWN
echo "On this machine the G4VIS_USE_DAWN=\$G4VIS_USE_DAWN"
fi 

#+
if [ X$g4vis_use_opacs = Xy ] ; then
G4VIS_USE_OPACS=1
export G4VIS_USE_OPACS
echo "On this machine the G4VIS_USE_OPACS=\$G4VIS_USE_OPACS"
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
if [ X$g4vis_oglhome != X ] ; then
OGLHOME="$g4vis_oglhome"
export OGLHOME
echo "On this machine the OGLHOME=\$OGLHOME"
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

if [ X$g4lib_use_shared = Xy ] ; then
G4LIB_USE_SHARED=1
export G4LIB_USE_SHARED
echo "On this machine the G4LIB_USE_SHARED=\$G4LIB_USE_SHARED"
fi 

if [ X$g4lib_use_static = Xy ] ; then
G4LIB_USE_STATIC=1
export G4LIB_USE_STATIC
echo "On this machine the G4LIB_USE_STATIC=\$G4LIB_USE_STATIC"
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

!GROK!THIS!

: In the following dollars and backticks do not need the extra backslash.
$spitshell >>env.sh <<'!NO!SUBS!'

#+
#
# G4WORKDIR
#
if [ $G4WORKDIR ] ; then
echo "In your environment you have the G4WORKDIR=$G4WORKDIR"
else
echo "G4WORKDIR will be set to HOME=$HOME."
export G4WORKDIR=$HOME
fi

!NO!SUBS!
chmod 755 env.sh
$eunicefix env.sh
