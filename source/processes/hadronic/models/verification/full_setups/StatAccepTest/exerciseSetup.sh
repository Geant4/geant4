#!/bin/sh
#
#----------------------------------------------------------------------------
# Last update: 19-Nov-2009.
#
# Setup script. 
#----------------------------------------------------------------------------
#
#***LOOKHERE***
export VO_GEANT4_SW_DIR=$PWD/..
#
export DIR_INSTALLATIONS=$VO_GEANT4_SW_DIR/dirInstallations
#
###export FULLY_STATIC_EXECUTABLE=1
#
# ===========
# === GCC ===
# ===========
#
###export PATH=$DIR_INSTALLATIONS/dirGCC/bin:$PATH
###export LD_LIBRARY_PATH=$DIR_INSTALLATIONS/dirGCC/lib:$LD_LIBRARY_PATH
#
###echo ' which g++ '
###which g++
###echo ' g++ -v '
###g++ -v
###echo ' '
#
# ==============
# === GEANT4 ===
# ==============
#
export G4SYSTEM=Linux-g++
#
export G4_RELEASE=dirGeant4-9.2.p02
#
export G4INSTALL=$DIR_INSTALLATIONS/$G4_RELEASE
export G4LIB=$G4INSTALL/lib
#
export G4LEVELGAMMADATA=$DIR_INSTALLATIONS/data/PhotonEvaporation
export G4RADIOACTIVEDATA=$DIR_INSTALLATIONS/data/RadioactiveDecay
export G4LEDATA=$DIR_INSTALLATIONS/data/G4EMLOW
export NeutronHPCrossSections=$DIR_INSTALLATIONS/data/G4NDL
export G4NEUTRONHPDATA=$DIR_INSTALLATIONS/data/G4NDL
#
export CLHEP_BASE_DIR=$DIR_INSTALLATIONS/dirCLHEP-2.0.4.2
export CLHEP_INCLUDE_DIR=$CLHEP_BASE_DIR/include
export CLHEP_LIB_DIR=$CLHEP_BASE_DIR/lib
export CLHEP_LIB=CLHEP
#
export G4UI_USE_TCSH=1
export G4VIS_NONE=1
#
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$G4LIB/$G4SYSTEM
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CLHEP_LIB_DIR
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$G4LIB/plists/$G4SYSTEM
#
# ===================
# === APPLICATION ===
# ===================
#
# --- Setup your Geant4 application environment ---
export G4WORKDIR=$PWD
export PATH=$PATH:$G4WORKDIR/bin/$G4SYSTEM
export G4ANALYSIS_USEROOT=1
#
# --- ROOT ---
export ROOTSYS=$DIR_INSTALLATIONS/dirROOT
export PATH=$PATH:$ROOTSYS/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTSYS/lib
#
export PHYSLIST=FTFP
