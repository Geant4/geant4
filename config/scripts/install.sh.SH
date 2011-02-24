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
# Clean all G4 envs
unset  CLHEP_BASE_DIR
unset  CLHEP_INCLUDE_DIR
unset  CLHEP_LIB
unset  CLHEP_LIB_DIR

unset  G4DEBUG
unset  G4INCLUDE
unset  G4INSTALL

unset  G4LEDATA
unset  G4LEVELGAMMADATA
unset  G4NEUTRONHPDATA
unset  G4RADIOACTIVEDATA
unset  G4ABLADATA

unset  G4LIB
unset  G4LIB_BUILD_G3TOG4
unset  G4LIB_BUILD_SHARED
unset  G4LIB_BUILD_STATIC
unset  G4LIB_BUILD_ZLIB
unset  G4LIB_BUILD_GDML
unset  G4LIB_USE_G3TOG4
unset  G4LIB_USE_GRANULAR
unset  G4LIB_USE_ZLIB

unset  G4SYSTEM

unset  G4UI_NONE
unset  G4UI_BUILD_WIN32_SESSION
unset  G4UI_BUILD_XAW_SESSION
unset  G4UI_BUILD_XM_SESSION
unset  G4UI_BUILD_QT_SESSION
unset  G4UI_USE_TCSH
unset  G4UI_USE_WIN32
unset  G4UI_USE_XAW
unset  G4UI_USE_XM
unset  G4UI_USE_QT

unset  G4VIS_NONE
unset  G4VIS_BUILD_DAWN_DRIVER
unset  G4VIS_BUILD_OIWIN32_DRIVER
unset  G4VIS_BUILD_OIX_DRIVER
unset  G4VIS_BUILD_OPENGLWIN32_DRIVER
unset  G4VIS_BUILD_OPENGLXM_DRIVER
unset  G4VIS_BUILD_OPENGLX_DRIVER
unset  G4VIS_BUILD_RAYTRACERX_DRIVER
unset  G4VIS_BUILD_VRML_DRIVER
unset  G4VIS_BUILD_OPENGLQT_DRIVER

unset  G4VIS_USE_DAWN
unset  G4VIS_USE_OIWIN32
unset  G4VIS_USE_OIX
unset  G4VIS_USE_OPENGLWIN32
unset  G4VIS_USE_OPENGLX
unset  G4VIS_USE_OPENGLXM
unset  G4VIS_USE_RAYTRACERX
unset  G4VIS_USE_VRML
unset  G4VIS_USE_OPENGLQT

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
#echo "On this machine the G4ELASTICDATA=\$g4elasticdata"
#fi

#
# g4clhep.U
#
if [ X$g4clhep_base_dir != X ] ; then 
CLHEP_BASE_DIR="$g4clhep_base_dir"
export CLHEP_BASE_DIR
echo "On this machine the CLHEP_BASE_DIR=\$CLHEP_BASE_DIR"
fi

#if [ X$g4clhep_base_dir != X/usr -o X$g4clhep_base_dir != X/usr/local ] ; then 
#CLHEP_BASE_DIR="$g4clhep_base_dir"
#export CLHEP_BASE_DIR
#echo "On this machine the CLHEP_BASE_DIR=\$CLHEP_BASE_DIR"
#fi

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
if [ X$g4ui_none = Xy ] ; then
G4UI_NONE=1
export G4UI_NONE
echo "On this machine G4UI_NONE is set - no UIs will be built"

else
echo "On this machine G4UI_NONE is not set - following UIs will be built:"

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

fi
# END OF G4UI_NONE block

#
# g4vis
#
#+
if [ X$g4vis_none = Xy ] ; then
G4VIS_NONE=1
export G4VIS_NONE
echo "On this machine G4VIS_NONE is set - no drivers will be built"

else
echo "On this machine G4VIS_NONE is not set - following drivers will be built:"

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
if [ X$g4vis_build_raytracerx_driver = Xy ] ; then
G4VIS_BUILD_RAYTRACERX_DRIVER=1
export G4VIS_BUILD_RAYTRACERX_DRIVER
echo "On this machine the G4VIS_BUILD_RAYTRACERX_DRIVER=\$G4VIS_BUILD_RAYTRACERX_DRIVER"
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
if [ X$g4vis_use_raytracerx = Xy ] ; then
G4VIS_USE_RAYTRACERX=1
export G4VIS_USE_RAYTRACERX
echo "On this machine the G4VIS_USE_RAYTRACERX=\$G4VIS_USE_RAYTRACERX"
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


fi
# End of G4VIS_NONE block



#+
if [ X$g4vis_oglhome != X ] ; then
OGLHOME="$g4vis_oglhome"
export OGLHOME
echo "On this machine the OGLHOME=\$OGLHOME"
fi 

#+
if [ X$g4vis_oivhome != X ] ; then
OIVHOME="$g4vis_oivhome"
export OIVHOME
if [ X\$g4non_display = X ] ; then
echo "On this machine the OIVHOME=\$OIVHOME"
fi
fi 


#+
if [ X$g4vis_build_openglxm_driver != X ] ; then
XMFLAGS="$g4vis_xmflags"
export XMFLAGS
if [ X\$g4non_display = X ] ; then
echo "On this machine the XMFLAGS=\$XMFLAGS"
fi
fi 

#+
if [ X$g4vis_build_openglxm_driver != X ] ; then
XMLIBS="$g4vis_xmlibs"
export XMLIBS
if [ X\$g4non_display = X ] ; then
echo "On this machine the XMLIBS=\$XMLIBS"
fi
fi 

#+
if [ X$g4ui_build_xm_session != X ] ; then
XMFLAGS="$g4vis_xmflags"
export XMFLAGS
if [ X\$g4non_display = X ] ; then
echo "On this machine the XMFLAGS=\$XMFLAGS"
fi
fi 

#+
if [ X$g4vis_build_xm_session != X ] ; then
XMLIBS="$g4vis_xmlibs"
export XMLIBS
if [ X\$g4non_display = X ] ; then
echo "On this machine the XMLIBS=\$XMLIBS"
fi
fi 

#+
if [ X$g4ui_build_xaw_session != X ] ; then
XAWFLAGS="$g4ui_xawflags"
export XAWFLAGS
if [ X\$g4non_display = X ] ; then
echo "On this machine the XAWFLAGS=\$XAWFLAGS"
fi
fi 

#+
if [ X$g4ui_build_xaw_session != X ] ; then
XAWLIBS="$g4ui_xawlibs"
export XAWLIBS
if [ X\$g4non_display = X ] ; then
echo "On this machine the XAWLIBS=\$XAWLIBS"
fi
fi 



# Qt Flags and Libs, messy, but needed for backward compatibility
#+
if [ "X$g4ui_build_qt_session" = "Xy" -o "X$g4ui_use_qt" = "Xy" ] ; then
QTFLAGS="$g4_qt_cxxflags"
QTLIBS="$g4_qt_libs"
QTMOC="$g4_qt_moc"
export QTFLAGS
export QTLIBS
export QTMOC
fi

if [ "X$g4vis_build_openglqt_driver" = "Xy" -o "X$g4vis_build_openglqt_driver" = "Xy" ] ; then
if [ "X\$QTFLAGS" = "X" ] ; then
QTFLAGS="$g4_qt_cxxflags $g4_qt_gl_cxxflags"
else
QTFLAGS="\$QTFLAGS $g4_qt_gl_cxxflags"
fi
if [ "X\$QTMOC" = "X" ] ; then
QTMOC="$g4_qt_moc"
fi
GLQTLIBS="$g4_qt_libs $g4_qt_gl_libs"
export QTFLAGS
export QTMOC
export GLQTLIBS
fi

#
# Because the setting of flags is messy, we separate out the reporting
#
if test "x\$QTMOC" != "x" -a "x\$g4non_display" = "x" ; then
    echo "On this machine the QTMOC=\$QTMOC"
    echo "On this machine the QTFLAGS=\$QTFLAGS"
    # {GL}QTLIBS dependent on what modules were chosen...
    if test "x\$QTLIBS" != "x" ; then
        echo "On this machine the QTLIBS=\$QTLIBS"
    fi
    if test "x\$GLQTLIBS" != "x" ; then
        echo "On this machine the GLQTLIBS=\$GLQTLIBS"
    fi
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
# Create GNUmakefile.db to keep ALL GNUmakefile DB :)

if test -f \$G4INSTALL/config/scripts/GNUmakefile.db ; then
    mv \$G4INSTALL/config/scripts/GNUmakefile.db \$G4LIB/\$G4SYSTEM/GNUmakefile.db.bac
fi

# Be careful: first time there is no $G4LIB!

mkdir -p \$G4LIB/\$G4SYSTEM

$g4make -n -p -f \$G4INSTALL/config/scripts/GNUmakefile > \$G4LIB/\$G4SYSTEM/GNUmakefile.db

echo ""
echo "Starting build..."
echo ""

echo "Cleaning up $g4conf/g4make.log"
rm -rf $g4conf/g4make.log

cd $g4install/source

if test $g4make_jobs -gt 1 ; then
    $g4make obj -j $g4make_jobs
fi

# There is no env!
if [ X$g4global = Xy ] ; then
    echo "Building global libraries..."
    $g4make global 2>&1 | tee $g4conf/g4make.log
fi

# There is env for binaries - for User: should be set in env.[c]sh!
if [ X$g4granular = Xy ] ; then
    echo "Building granular libraries..."
    G4LIB_USE_GRANULAR=1
    export G4LIB_USE_GRANULAR
    $g4make 2>&1 | tee $g4conf/g4make.log
fi

if [ X$g4lib_build_dll = Xy ] ; then
    echo "Building WIN32 DLLs..."
    $g4make dll 2>&1 | tee $g4conf/g4make.log
fi

if [ X$g4lib_build_dyn = Xy ] ; then
    echo ""
    rm -rf $g4conf/g4make.log
    $g4make global G4LIB_BUILD_SHARED=1 2>&1 | tee $g4conf/g4make.log
fi

# Check for errors in log file
if test -f \$g4conf/g4make.log ; then
    # There are two (at least) shells 'sh' and 'bash' with
    # different arithmetic expansion...So we are doing it by 'text'
    # in two steps...
 
    err=\`grep " Error " $g4conf/g4make.log | wc -w\` 
    int=\`grep " Interrupt " $g4conf/g4make.log | wc -w\` 
    stp=\`grep Stop\. $g4conf/g4make.log | wc -w\` 

    err2=\`echo \$err\`
    int2=\`echo \$int\`
    stp2=\`echo \$stp\`

    if [ X\$err2 = X0 -a X\$int2 = X0 ] ; then
        echo ""
        echo "####################################################"
        echo "# Your Geant4 installation seems to be successful!  "
        echo "# To be sure please have a look into the log file:  "
        echo "# $g4conf/g4make.log				  "
        echo "####################################################"
        echo ""
    else
        echo ""
        echo "############################################################"
        echo "# It looks like you had errors during Geant4 installation   "
        echo "# (or interrupted installation)!				  "
        echo "# Please check log file:				      	  "
        echo "# $g4conf/g4make.log					  "
        echo "############################################################"
        echo ""
    fi 
else
    echo ""
    echo "############################################################"
    echo "# It looks like you had errors during Geant4 installation   "
    echo "# (or interrupted installation)!				  "
    echo "# You have no log file:				      	  "
    echo "# $g4conf/g4make.log					  "
    echo "############################################################"
    echo ""
fi

!GROK!THIS!

: In the following dollars and backticks do not need the extra backslash.
$spitshell >>install.sh <<'!NO!SUBS!'
!NO!SUBS!
chmod 755 install.sh
$eunicefix install.sh

