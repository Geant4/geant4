setenv G4INSTALL $HOME/geant4/cvsco/geant4
setenv G4SYSTEM Linux-g++
setenv G4LIB_BUILD_SHARED 1

# CLHEP (1.80)
set clhep_vers = 1.8.0.0
setenv CLHEP_BASE_DIR /usr/local/CLHEP/$clhep_vers

# ROOT
setenv ROOTSYS $HOME/local/root/v3.03.06/root
setenv PATH $ROOTSYS/bin:$PATH
setenv LD_LIBRARY_PATH $ROOTSYS/lib:$LD_LIBRARY_PATH

