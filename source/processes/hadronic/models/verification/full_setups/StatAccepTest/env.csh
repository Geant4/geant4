#
# ================
# === gcc 3.23 ===
# ================
#
source /afs/cern.ch/sw/geant4/dev/scripts/gcc-alt.csh 3.2.3
#
# ==============
# === GEANT4 ===
# ==============
#
setenv G4SYSTEM Linux-g++
#
setenv RELEASE geant4.6.2.ref01
setenv PLATFORM rh73_gcc323/
setenv DIR_SPECIFIC /afs/cern.ch/sw/geant4/releases/specific/$PLATFORM
setenv G4INSTALL /afs/cern.ch/sw/geant4/releases/share/$RELEASE
setenv G4LIB $DIR_SPECIFIC/$RELEASE/lib
#
setenv G4LEVELGAMMADATA /afs/cern.ch/sw/geant4/dev/data/PhotonEvaporation
setenv G4RADIOACTIVEDATA /afs/cern.ch/sw/geant4/dev/data/RadiativeDecay
setenv G4LEDATA /afs/cern.ch/sw/geant4/dev/data/G4EMLOW
setenv NeutronHPCrossSections /afs/cern.ch/sw/geant4/dev/data/G4NDL
#
setenv CLHEP_BASE_DIR /afs/cern.ch/sw/geant4/dev/CLHEP/1.8.1.0/redhat73_gcc323
setenv CLHEP_INCLUDE_DIR $CLHEP_BASE_DIR/include
setenv CLHEP_LIB_DIR $CLHEP_BASE_DIR/lib
setenv CLHEP_LIB CLHEP
#
setenv G4ANALYSIS_USE 1
setenv G4UI_USE_TCSH 1
#
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:$G4LIB/$G4SYSTEM
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:$G4LIB/.lists_build/$G4SYSTEM
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:$CLHEP_LIB_DIR
#
# ===================
# === APPLICATION ===
# ===================
#
# --- Setup your Geant4 application environment ---
setenv G4WORKDIR $PWD
setenv PATH ${PATH}:$G4WORKDIR/bin/$G4SYSTEM
setenv G4ANALYSIS_USE 1
#
# --- Aida / PI ---
setenv PI_DIR /afs/cern.ch/sw/lcg/app/releases/PI/PI_1_2_3/rh73_gcc323
setenv PATH $PI_DIR/bin:${PATH}
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:$PI_DIR/lib
eval `aida-config --runtime csh`
#
#***LOOKHERE*** TEMPORARY LORENZO FIX: CLOUD CACHE SETS TO 10,000 .
setenv LD_LIBRARY_PATH /afs/cern.ch/user/r/ribon/ExtraSpace/myPiDir/lib:${LD_LIBRARY_PATH}
#
# --- GSL : this is needed only for  dirStat/pvalue.cpp ---
setenv GSL_DIR /afs/cern.ch/sw/lcg/external/GSL/1.4/rh73_gcc323
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:$GSL_DIR/lib
#

