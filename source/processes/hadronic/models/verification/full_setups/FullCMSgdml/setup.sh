#!/bin/sh
#
#-----------------------------------------------------------------
# Last update: 28-Nov-2007
#-----------------------------------------------------------------
#
# ==============================
# === Special setup for GDML ===
# ==============================
#
export G4LIB_BUILD_GDML=1
export G4LIB_USE_GDML=1
export XERCESCROOT=/afs/cern.ch/sw/geant4/dev/XercesC/2.8.0/slc4_gcc346
export LD_LIBRARY_PATH=$XERCESCROOT/lib:$LD_LIBRARY_PATH
#
# ===========================================
# === Usual setup for Geant4 applications ===
# ===========================================
#
#***LOOKHERE*** Choose the Geant4 and CLHEP versions
#
# --- G4 6.2.p02 ---
###g4version="6.2.p02"
###clhepversion="1.8.1.0"
#
# --- G4 7.1.p01a ---
###g4version="7.1.p01a"
###clhepversion="1.9.2.2"
#
# --- G4 8.0.p01 ---
###g4version="8.0.p01"
###clhepversion="1.9.2.2"
#
# --- G4 8.1.p02 ---
###g4version="8.1.p02"
###clhepversion="1.9.2.3"
#
# --- G4 8.2.p01 ---
###g4version="8.2.p01"
###clhepversion="2.0.3.1"
#
# --- G4 9.0.p01 --
g4version="geant4.9.0.p01"
clhepversion="2.0.3.1"
#
#**************
#
g4releases=/afs/cern.ch/sw/lcg/external/geant4
clhep=/afs/cern.ch/sw/lcg/external/clhep
#
gccversion="3.2.3"
os="slc3_ia32_gcc323"
#
g++ --version | grep $gccversion > /dev/null
if [ $? != 0 ]
then
  echo "It looks like your compiler settings are not suitable"
  echo "The Operating system is expected to be $os"
  echo    "The compiler version should be g++ (GCC) $gccversion"
  echo -n "The system reports that it is  "; g++ --version
  echo "Please set your PATH and LD_LIBRARY_PATH environment variables"
#
else
#
echo "Setting up the environment for $g4version"
#
G4SYSTEM=Linux-g++
G4INSTALL=$g4releases/$g4version/share
G4LIB=$g4releases/$g4version/$os/lib
CLHEP_BASE_DIR=$clhep/$clhepversion/$os
export G4SYSTEM G4INSTALL G4LIB CLHEP_BASE_DIR
#
G4LEDATA="$G4INSTALL/data/G4EMLOW"
NeutronHPCrossSections="$G4INSTALL/data/G4NDL"
G4LEVELGAMMADATA="$G4INSTALL/data/PhotonEvaporation"
G4RADIOACTIVEDATA="$G4INSTALL/data/RadiativeDecay"
G4ELASTICDATA="$G4INSTALL/data/G4ELASTIC"
export G4LEDATA NeutronHPCrossSections G4LEVELGAMMADATA G4RADIOACTIVEDATA G4ELASTICDATA
#
# Geant 4 interface, visualisation and other variables
G4UI_USE_TERMINAL=1
G4UI_USE_TCSH=1
G4UI_USE_GAG=1
G4UI_USE_XAW=1
G4UI_USE_XM=1
#
G4VIS_USE_DAWN=1
G4VIS_USE_DAWNFILE=1
G4VIS_USE_OPENGLX=1
G4VIS_USE_OPENGLXM=1
G4VIS_USE_RAYTRACER=1
G4VIS_USE_RAYTRACERX=1
G4VIS_USE_VRML=1
G4VIS_USE_VRMLFILE=1
#
G4LIB_USE_G3TOG4=1
#
export G4UI_USE_TERMINAL G4UI_USE_TCSH G4UI_USE_GAG G4UI_USE_XAW G4UI_USE_XM
export G4VIS_USE_DAWN G4VIS_USE_DAWNFILE G4VIS_USE_OPENGLX G4VIS_USE_OPENGLXM
export G4VIS_USE_RAYTRACER G4VIS_USE_RAYTRACERX G4VIS_USE_VRML G4VIS_USE_VRMLFILE
export G4LIB_USE_G3TOG4
#
if [ $LD_LIBRARY_PATH ]
then
    LD_LIBRARY_PATH=${G4LIB}:${G4LIB}/plists/${G4SYSTEM}:${CLHEP_BASE_DIR}/lib:${LD_LIBRARY_PATH}
else
    LD_LIBRARY_PATH=${G4LIB}:${G4LIB}/plists/${G4SYSTEM}:${CLHEP_BASE_DIR}/lib
fi
export LD_LIBRARY_PATH
#
fi
#
# ===================
# === Application ===
# ===================
#
export G4WORKDIR=$PWD
export G4TMP=$G4WORKDIR/tmp
export G4BIN=$G4WORKDIR/bin
#
export PATH=$PATH:$G4WORKDIR/bin/$G4SYSTEM
#
###export G4ANALYSIS_USE=1
#
