#!/usr/local/bin/bash
export G4INSTALL=/users/sscherer/geant4
cd $G4INSTALL/source
pwd
echo ">>>> gmake clean; gmake"
gmake clean; gmake
cd $G4INSTALL/source/processes/hadronic/models/quark_molecular_dynamics/util/
pwd
echo ">>>> gmake clean; gmake"
gmake clean; gmake
cd $G4INSTALL/source/processes/hadronic/models/quark_molecular_dynamics/body/
pwd
echo ">>>> gmake clean; gmake"
gmake clean; gmake
cd $G4INSTALL/source
pwd
echo ">>>> gmake libmap"
gmake libmap
cd $G4INSTALL/source/processes/hadronic/models/quark_molecular_dynamics/test/
pwd
export G4TARGET=G4qmdTest
echo ">>>> gmake with G4TARGET=$G4TARGET"
gmake
