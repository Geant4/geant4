#!/bin/sh
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
export G4_RELEASE=dirGeant4-9.2.p01
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
export CLHEP_BASE_DIR=$DIR_INSTALLATIONS/dirCLHEP
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
export G4ANALYSIS_USE=1
#
# --- AIDA ---
export CERNLIB_DIR=$DIR_INSTALLATIONS/dirCERNLIB
export PATH=$DIR_INSTALLATIONS/diriAIDA/bin:$PATH
eval `aida-config --runtime sh`
#
# --- Statistical toolkit ---
export DIR_STAT=$DIR_INSTALLATIONS/dirStatisticalToolkit
export LD_LIBRARY_PATH=$DIR_STAT/lib:$LD_LIBRARY_PATH
#
# --- GSL : this is needed only for  dirStat/pvalue.cpp ---
export GSL_DIR=$DIR_INSTALLATIONS/dirGSL
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GSL_DIR/lib
#
# --- PAW
export PATH=$PATH:$DIR_INSTALLATIONS/dirPAW
#
export PHYSLIST=FTFP
