#
# ================
# === gcc 3.23 ===
# ================
#
. /afs/cern.ch/sw/geant4/dev/scripts/gcc-alt.sh 3.2.3
#
# ==============
# === GEANT4 ===
# ==============
#
export G4SYSTEM=Linux-g++
#
RELEASE=geant4.6.2.ref03
PLATFORM=rh73_gcc323/
DIR_SPECIFIC=/afs/cern.ch/sw/geant4/releases/specific/$PLATFORM
export G4INSTALL=/afs/cern.ch/sw/geant4/releases/share/$RELEASE
export G4LIB=$DIR_SPECIFIC/$RELEASE/lib
#
export G4LEVELGAMMADATA=/afs/cern.ch/sw/geant4/dev/data/PhotonEvaporation
export G4RADIOACTIVEDATA=/afs/cern.ch/sw/geant4/dev/data/RadiativeDecay
export G4LEDATA=/afs/cern.ch/sw/geant4/dev/data/G4EMLOW
export NeutronHPCrossSections=/afs/cern.ch/sw/geant4/dev/data/G4NDL
#
export CLHEP_BASE_DIR=/afs/cern.ch/sw/geant4/dev/CLHEP/1.8.1.0/redhat73_gcc323
export CLHEP_INCLUDE_DIR=$CLHEP_BASE_DIR/include
export CLHEP_LIB_DIR=$CLHEP_BASE_DIR/lib
export CLHEP_LIB=CLHEP
#
export G4ANALYSIS_USE=1
export G4UI_USE_TCSH=1
#
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$G4LIB/$G4SYSTEM
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$G4LIB/.lists_build/$G4SYSTEM
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CLHEP_LIB_DIR
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
# --- Aida / PI ---
export PI_DIR=/afs/cern.ch/sw/lcg/app/releases/PI/PI_1_2_4/rh73_gcc323
export PATH=$PI_DIR/bin:$PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PI_DIR/lib
eval `aida-config --runtime sh`
#
# --- GSL : this is needed only for  dirStat/pvalue.cpp ---
export GSL_DIR=/afs/cern.ch/sw/lcg/external/GSL/1.4/rh73_gcc323
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GSL_DIR/lib
#
