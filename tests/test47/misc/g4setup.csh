#!/bin/tcsh -f

# setup all G4 variables manually
#
setenv G4SYSTEM Linux-g++
setenv G4INSTALL /uscmst1/prod/sw/cms/slc4_ia32_gcc345/external/geant4/9.2
setenv G4SRC     $G4INSTALL/source
setenv G4LIB     $G4INSTALL/lib
setenv G4INCLUDE $G4INSTALL/include

setenv G4DATA $G4INSTALL/data

# G4 service stuff
#
setenv G4LIB_BUILD_SHARED 1
setenv GLOBALLIBS 1

setenv CLHEP_BASE_DIR /uscmst1/prod/sw/cms/slc4_ia32_gcc345/external/clhep/1.9.4.2
setenv CLHEP_LIB CLHEP

setenv G4LEVELGAMMADATA  $G4DATA/PhotonEvaporation2.0
setenv G4NEUTRONHPDATA   $G4DATA/G4NDL3.13
setenv G4RADIOACTIVEDATA $G4DATA/RadioactiveDecay3.2
setenv G4LEDATA          $G4DATA/G4EMLOW6.5
setenv G4ABLADATA        $G4DATA/G4ABLA3.0


setenv ROOTSYS /uscmst1/prod/sw/cms/slc4_ia32_gcc345/lcg/root/5.22.00a-cms4/

setenv G4WORKDIR  $PWD

setenv G4EXE $G4WORKDIR/bin/$G4SYSTEM

setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH\:$G4LIB/:$CLHEP_BASE_DIR/lib/:$ROOTSYS/lib/
