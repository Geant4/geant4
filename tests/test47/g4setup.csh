#!/usr/local/bin/tcsh

setup ups
setup upd

ups list -aK+ geant4 v4_9_2
setup geant4 v4_9_2  -q g77-OpenGL 

ups list -aK+ root   v5_22_00
setup root v5_18_00  -q GCC_3_4_3

setenv G4WORKDIR  $PWD

setenv G4EXE $G4WORKDIR/bin/$G4SYSTEM

setup g4neutron   v3_13 
setenv G4NEUTRONHPDATA $NeutronHPCrossSections 

setup g4photon    v2_0
  
setup g4emlow     v6_2
    
setup g4radiative v3_2
    
setup g4elastic   v1_1

# this is optional; a typical G4 application can run without it 
#
###setup g4abla V3_0

setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH\:$G4LIB/plists/Linux-g++/:$G4LIB/Linux-g++/:$CLHEP_BASE_DIR/lib

