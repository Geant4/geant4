#!/bin/sh -l
set -x
export TOP=$(pwd)
brew tap davidchall/hep
brew install hepmc lhapdf wget cmake coreutils  
brew install gnu-sed
brew install gcc
brew install expat
brew install pythia
brew install libxmu
brew install --cask xquartz

which gfortran-11
if [ "$?" = "0" ]; then 
   export CXX=g++-11
   export CC=gcc-11
   export FC=gfortran-11
   export F77=gfortran-11
   export LD=gfortran-11
else
   export CXX=g++
   export CC=gcc
   export FC=gfortran
   export F77=gfortran
   export LD=gfortran
fi
export CXX=clang++
export CC=clang
cp /usr/local/bin/gfortran-11 /usr/local/bin/gfortran

ls -lah
wget -q  https://geant4-data.web.cern.ch/releases/lib_11.0.1/Darwin-clang13.0.0-Monterey.tar.gz
tar -xzf Darwin-clang13.0.0-Monterey.tar.gz
ls -lah

mkdir -p $TOP/Geant4-11.0.1-Darwin/share/Geant4-11.0.1/data/
cd $TOP/Geant4-11.0.1-Darwin/share/Geant4-11.0.1/data/
wget -q  https://geant4-data.web.cern.ch/datasets/G4NDL.4.6.tar.gz
wget -q  https://geant4-data.web.cern.ch/datasets/G4EMLOW.8.0.tar.gz
wget -q  https://geant4-data.web.cern.ch/datasets/G4PhotonEvaporation.5.7.tar.gz
wget -q  https://geant4-data.web.cern.ch/datasets/G4RadioactiveDecay.5.6.tar.gz
wget -q  https://geant4-data.web.cern.ch/datasets/G4PARTICLEXS.4.0.tar.gz
wget -q  https://geant4-data.web.cern.ch/datasets/G4PII.1.3.tar.gz
wget -q  https://geant4-data.web.cern.ch/datasets/G4RealSurface.2.2.tar.gz
wget -q  https://geant4-data.web.cern.ch/datasets/G4SAIDDATA.2.0.tar.gz
wget -q  https://geant4-data.web.cern.ch/datasets/G4ABLA.3.1.tar.gz
wget -q  https://geant4-data.web.cern.ch/datasets/G4INCL.1.0.tar.gz
wget -q  https://geant4-data.web.cern.ch/datasets/G4ENSDFSTATE.2.3.tar.gz
 
tar -xzf G4NDL.4.6.tar.gz
tar -xzf G4EMLOW.8.0.tar.gz
tar -xzf G4PhotonEvaporation.5.7.tar.gz
tar -xzf G4RadioactiveDecay.5.6.tar.gz
tar -xzf G4PARTICLEXS.4.0.tar.gz
tar -xzf G4PII.1.3.tar.gz
tar -xzf G4RealSurface.2.2.tar.gz
tar -xzf G4SAIDDATA.2.0.tar.gz
tar -xzf G4ABLA.3.1.tar.gz
tar -xzf G4INCL.1.0.tar.gz
tar -xzf G4ENSDFSTATE.2.3.tar.gz

ls -lah
cd $TOP


gsed -i 's@/Users/gcosmo/Software/release/install/@'$TOP'/Geant4-11.0.1-Darwin/@g'  $TOP/Geant4-11.0.1-Darwin/bin/geant4.sh

source $TOP/Geant4-11.0.1-Darwin/bin/geant4.sh
gsed -i 's@/opt/local/@/usr/local/Cellar/expat/2.4.7/@g'  $TOP/Geant4-11.0.1-Darwin/lib/Geant4-11.0.1/Geant4PackageCache.cmake

mkdir -p $TOP/test
cmake -B $TOP/test -S $TOP/examples/extended/eventgenerator -DPYTHIA6_INTERNAL=ON -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON -DCMAKE_INSTALL_PREFIX=$TOP/INSTALL -DGeant4_DIR=$TOP/Geant4-11.0.1-Darwin/lib/Geant4-11.0.1  -DHEPMC_DIR=/usr/local -DCMAKE_Fortran_FLAGS=-fPIC -DCMAKE_C_FLAGS=-fPIC -DCMAKE_CXX_FLAGS=-fPIC
cmake --build $TOP/test -j 2
cmake --install $TOP/test


#Running tests
pwd
find $TOP | grep hepmc_pygen.in
cd  $TOP/examples/extended/eventgenerator/HepMC/HepMCEx01
$TOP/INSTALL/bin/HepMCEx01 hepmc_pygen.in
$TOP/INSTALL/bin/HepMCEx01 hepmc_ascii.in

cd  $TOP/examples/extended/eventgenerator/HepMC/HepMCEx02
$TOP/INSTALL/bin/HepMCEx02 hepmc_pygen.in
$TOP/INSTALL/bin/HepMCEx02 hepmc_ascii.in
