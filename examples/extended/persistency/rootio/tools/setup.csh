setenv G4INSTALL $HOME/geant4/cvsco/geant4
setenv G4SYSTEM Linux-g++
setenv G4LIB_BUILD_SHARED 1

# CLHEP (1.80)
set clhep_vers = 1.8.0.0
setenv CLHEP_BASE_DIR /usr/local/CLHEP/$clhep_vers

# ROOT
setenv ROOTSYS $HOME/KEK/local/root/3.03.08
setenv PATH $ROOTSYS/bin:$PATH

