#!/usr/local/bin/tcsh

setup ups
setup upd

ups list -aK+ geant4 | grep g4_9
###setup geant4 v4_9_2_p01  -q g77-OpenGL 
setup geant4 v4_9_4_b01

ups list -aK+ root   
setup root v5_18_00f  -q GCC_4_1_2

setenv G4WORKDIR  $PWD

setenv G4EXE $G4WORKDIR/bin/$G4SYSTEM

setup g4neutron   v3_13 
###setenv NeutronHPCrossSections g4neutron
setenv G4NEUTRONHPDATA $NeutronHPCrossSections 

setup g4photon    v2_0
###setenv G4LEVELGAMMADATA g4photon 
  
setup g4emlow     v6_9
###setenv G4LEDATA g4emlow
    
setup g4radiative v3_2
###setenv G4RADIOACTIVEDATA g4radiactive
    
setup g4elastic   v1_1
###setenv G4ELASTICDATA g4elastic 

setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH\:$G4LIB/$G4SYSTEM/:$CLHEP_BASE_DIR/lib/:$ROOTSYS/lib
