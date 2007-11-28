#
# ==============
# === GEANT4 ===
# ==============
#
export G4SYSTEM=Linux-g++
#
# Use of my own geant4 installation 
export G4INSTALL=/users/ribon/geant4
export G4LIB=$G4INSTALL/lib
#
# For official Geant4 releases, use  /afs/cern.ch/sw/lcg/external/geant4
# For developer Geant4 reference tags, use /afs/cern.ch/sw/geant4/releases
#
###export G4INSTALL=/afs/cern.ch/sw/lcg/external/geant4/7.1/share
###export G4LIB=/afs/cern.ch/sw/lcg/external/geant4/7.1/slc3_ia32_gcc323/lib
#
###export G4INSTALL=/afs/cern.ch/sw/geant4/releases/share/geant4.7.1.ref03
###export G4LIB=/afs/cern.ch/sw/geant4/releases/specific/slc3_ia32_gcc323/geant4.7.1.ref03/lib
#
export G4DEV=/afs/cern.ch/sw/geant4/dev
###export G4DEBUG=1
#
#--- data ---
export G4LEVELGAMMADATA=$G4DEV/data/PhotonEvaporation
export G4RADIOACTIVEDATA=$G4DEV/data/RadiativeDecay
export G4LEDATA=$G4DEV/data/G4EMLOW
export NeutronHPCrossSections=$G4DEV/data/G4NDL
#
#--- CLHEP ---
# Before Geant4 version 7.1, use CLHEP 1.9.1.2 ;
# for Geant4 version 7.1, use CLHEP 1.9.2.1 . 
# for Geant4 version 8.0, use CLHEP 1.9.2.2 .
# for Geant4 version 8.1, use CLHEP 1.9.2.3 .
###export CLHEP_BASE_DIR=$G4DEV/CLHEP/1.9.1.2/slc3_gcc323
###export CLHEP_BASE_DIR=$G4DEV/CLHEP/1.9.2.1/slc3_gcc323
###export CLHEP_BASE_DIR=$G4DEV/CLHEP/1.9.2.2/slc3_gcc323
export CLHEP_BASE_DIR=$G4DEV/CLHEP/1.9.2.3/slc3_gcc323
#
export CLHEP_INCLUDE_DIR=$CLHEP_BASE_DIR/include
export CLHEP_LIB_DIR=$CLHEP_BASE_DIR/lib
export CLHEP_LIB=CLHEP
#
# --- PATH ---
export PATH=$PATH:$G4BIN/$G4SYSTEM
#
# --- LD_LIBRARY_PATH ---
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$G4LIB/$G4SYSTEM
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CLHEP_LIB_DIR
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$G4LIB/plists/$G4SYSTEM
#
#
# ===================
# === APPLICATION ===
# ===================
#
# --- Setup your Geant4 application environment ---
#
export G4WORKDIR=$PWD
export G4TMP=$G4WORKDIR/tmp
export G4BIN=$G4WORKDIR/bin
#
export PATH=$PATH:$G4WORKDIR/bin/$G4SYSTEM
#
export G4VIS_NONE=1
#
export G4ANALYSIS_USE=1
#
# --- AIDA/PI ---
# Use PI 1.3.3 before G4 8.0; use 1.3.12 for G4 8.0 .
###export PI_DIR=$G4DEV/PI/1.3.3-lite/slc3_gcc323
export PI_DIR=$G4DEV/PI/1.3.12-lite/slc3_gcc323
export PATH=$PI_DIR/bin:$PATH
eval `aida-config --runtime sh`

