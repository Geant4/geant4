#!/usr/local/bin/tcsh

setup ups
setup upd

ups list -aK+ geant4 | grep g4_9
###setup geant4 v4_9_2_p01  -q g77-OpenGL 
setup geant4 v4_9_4_p01

ups list -aK+ root   
setup root v5_18_00f  -q GCC_4_1_2


setup g4neutron   v3_14 
###setenv NeutronHPCrossSections g4neutron
setenv G4NEUTRONHPDATA $NeutronHPCrossSections 

setup g4photon    v2_2
###setenv G4LEVELGAMMADATA g4photon 
  
setup g4emlow     v6_22
###setenv G4LEDATA g4emlow
    
setup g4radiative v3_3
###setenv G4RADIOACTIVEDATA g4radiactive
    
setup g4elastic   v1_1
###setenv G4ELASTICDATA g4elastic 


# example alternative setup
#
#setenv G4SYSTEM Linux-g++
#setenv G4INSTALL /uscmst1/prod/sw/cms/slc4_ia32_gcc345/external/geant4/9.2
#setenv G4SRC     $G4INSTALL/source
#setenv G4LIB     $G4INSTALL/lib
#setenv G4INCLUDE $G4INSTALL/include

# G4 service stuff
#
#setenv G4LIB_BUILD_SHARED 1
#setenv GLOBALLIBS 1

# may or may not be needed - as of g4.9.5, built-in clhep available
#
###setenv CLHEP_BASE_DIR /uscmst1/prod/sw/cms/slc4_ia32_gcc345/external/clhep/1.9.4.2
###setenv CLHEP_LIB CLHEP

#setenv G4DATA $G4INSTALL/data
#setenv G4LEVELGAMMADATA  $G4DATA/PhotonEvaporation2.2
#setenv G4NEUTRONHPDATA   $G4DATA/G4NDL4.0
#setenv G4RADIOACTIVEDATA $G4DATA/RadioactiveDecay3.3
#setenv G4LEDATA          $G4DATA/G4EMLOW6.22
#setenv G4ABLADATA        $G4DATA/G4ABLA3.0

#setenv ROOTSYS /uscmst1/prod/sw/cms/slc4_ia32_gcc345/lcg/root/5.22.00a-cms4/

# note: for cluck at fnal, poit to "central" installs 
# export G4DATA=/mnt/disk1/products
# export G4LEVELGAMMADATA=$G4DATA/g4photon/v2_2/NULL/PhotonEvaporation2.2
# export G4NEUTRONHPDATA=$G4DATA/g4neutron/v4_0/NULL/G4NDL4.0
# export G4RADIOACTIVEDATA=$G4DATA/g4radiative/v3_4/NULL/RadioactiveDecay3.4
# export G4LEDATA=$G4DATA/g4emlow/v6_23/NULL/G4EMLOW6.23
#
# export ROOTSYS=/mnt/disk1/products/root/v5_30_00/Linux64bit+2.6-2.5-gcc45


setenv G4WORKDIR  $PWD
setenv G4EXE $G4WORKDIR/bin/$G4SYSTEM

setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH\:$G4LIB/$G4SYSTEM/:$CLHEP_BASE_DIR/lib/:$ROOTSYS/lib
#
# for global libs, it's just $G4LIB path...
# and no need for clhep
#
#setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH\:$G4LIB/:$ROOTSYS/lib
#
# additional settings for testing nuclear structure parameters
#
# setenv G4NUCMODEL_XSEC_SCALE 0.1
# setenv G4NUCMODEL_RAD_SCALE 1.0
# setenv G4NUCMODEL_RAD_2PAR 1
# setenv G4NUCMODEL_RAD_SMALL 1.992
# setenv G4NUCMODEL_RAD_ALPHA 0.84
# setenv G4NUCMODEL_FERMI_SCALE 0.685
# setenv G4NUCMODEL_RAD_TRAILING 1.2
