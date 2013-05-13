#!/usr/bin/env bash
#
# export G4INSTALL=/home-cluck/syjun/g4p/work/g4.9.6.b01/geant4.9.6.b01
export G4INSTALL=/home/yarba_j/work/test-geant4.9.6-ref03/geant4-09-06-ref-03
export G4SYSTEM=Linux-g++
export G4SRC=${G4INSTALL}/source
export G4LIB=${G4INSTALL}/lib64
export G4INCLUDE=${G4INSTALL}/include/Geant4

export G4LIB_BUILD_SHARED=1
export GLOBALLIBS=1

# export G4DATA=${G4INSTALL}/data
export G4DATA=/home/g4p/pbs/download/g4data
export G4LEDATA=${G4DATA}/G4EMLOW6.32
export G4LEVELGAMMADATA=${G4DATA}/PhotonEvaporation2.3
export G4NEUTRONHPDATA=${G4DATA}/G4NDL4.3
export G4RADIOACTIVEDATA=${G4DATA}/RadioactiveDecay3.6
export G4ABLADATA=${G4DATA}/G4ABLA3.0
export G4REALSURFACEDATA=${G4DATA}/RealSurface1.0
export G4NEUTRONXSDATA=${G4DATA}/G4NEUTRONXS1.2
export G4PIIDATA=${G4DATA}/G4PII1.3
# export G4SAIDDATA=${G4DATA}/G4SAIDDATA1.1
export G4SAIDXSDATA=${G4DATA}/G4SAIDDATA1.1

# export G4WORKDIR=/home-cluck/yarba_j/work/test-central-geant4.9.6.b01/tests/test48-build
export G4WORKDIR=${G4INSTALL}/tests
export G4EXE=${G4WORKDIR}/bin/${G4SYSTEM}

#export ROOTSYS=/products/root/v5_30_02/Linux64bit+2.6-2.5-gcc46-prof
export ROOTSYS=/products/root/v5_30_06/Linux64bit+2.6-2.12-gcc46-prof

# export XERCESCROOT=/products/xerces_c/v3_1_1/Linux64bit+2.6-2.5-gcc46-prof

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${G4LIB}:${ROOTSYS}/lib   #:${XERCESCROOT}/lib

export PATH=${PATH}:${G4EXE}:${ROOTSYS}/bin

