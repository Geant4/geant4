#
# ================
# === gcc 3.23 ===
# ================
#
export PATH=$VO_DTEAM_SW_DIR/dirInstallations/dirGCC/bin:$PATH
export LD_LIBRARY_PATH=$VO_DTEAM_SW_DIR/dirInstallations/dirGCC/lib:$LD_LIBRARY_PATH
#
# ==============
# === GEANT4 ===
# ==============
#
export G4SYSTEM=Linux-g++
#
###export G4_RELEASE=geant4.6.2.p02
export G4_RELEASE=geant4.7.0.cand01
#
export DIR_INSTALLATIONS=$VO_DTEAM_SW_DIR/dirInstallations
export G4INSTALL=$DIR_INSTALLATIONS/$G4_RELEASE
export G4LIB=$DIR_INSTALLATIONS/$G4_RELEASE/lib
#
export G4LEVELGAMMADATA=$DIR_INSTALLATIONS/dirG4DATA/PhotonEvaporation
export G4RADIOACTIVEDATA=$DIR_INSTALLATIONS/dirG4DATA/RadiativeDecay
export G4LEDATA=$DIR_INSTALLATIONS/dirG4DATA/G4EMLOW
export NeutronHPCrossSections=$DIR_INSTALLATIONS/dirG4DATA/G4NDL
#
export CLHEP_BASE_DIR=$DIR_INSTALLATIONS/dirCLHEP
export CLHEP_INCLUDE_DIR=$CLHEP_BASE_DIR/include
export CLHEP_LIB_DIR=$CLHEP_BASE_DIR/lib
export CLHEP_LIB=CLHEP
#
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
# Use my own PI libraries instead of the AFS ones.
export PI_DIR=$DIR_INSTALLATIONS/dirPI
export LD_LIBRARY_PATH=$PI_DIR/lib:${LD_LIBRARY_PATH}
export SEAL_PLUGINS=$PI_DIR/lib/modules
#
# --- GSL : this is needed only for  dirStat/pvalue.cpp ---
### export GSL_DIR=$DIR_INSTALLATIONS/dirGSL
### export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GSL_DIR/lib
#
