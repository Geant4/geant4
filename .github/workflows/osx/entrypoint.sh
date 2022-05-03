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

mkdir -p pythia6428-split
cd pythia6428-split
wget https://pythia.org/download/pythia6/pythia6428-split.tgz
tar -xzf pythia6428-split.tgz
gsed -i 's/pdfset.f//g' Makefile
gsed -i 's/gfortran/gfortran-11/g' Makefile
make lib
mv libpythia.a /usr/local/lib/libpythia6.a
cd $TOP

ls -lah
wget https://geant4-data.web.cern.ch/releases/lib_11.0.1/Darwin-clang13.0.0-Monterey.tar.gz
tar -xzf Darwin-clang13.0.0-Monterey.tar.gz
ls -lah

mkdir -p Geant4-11.0.1-Darwin/share/Geant4-11.0.1/data/
cd Geant4-11.0.1-Darwin/share/Geant4-11.0.1/data/
wget https://geant4-data.web.cern.ch/datasets/G4NDL4.6.tar.gz
wget https://geant4-data.web.cern.ch/datasets/G4EMLOW8.0.tar.gz
wget https://geant4-data.web.cern.ch/datasets/PhotonEvaporation5.7.tar.gz
wget https://geant4-data.web.cern.ch/datasets/RadioactiveDecay5.6.tar.gz
wget https://geant4-data.web.cern.ch/datasets/G4PARTICLEXS4.0.tar.gz
wget https://geant4-data.web.cern.ch/datasets/G4PII1.3.tar.gz
wget https://geant4-data.web.cern.ch/datasets/RealSurface2.2.tar.gz
wget https://geant4-data.web.cern.ch/datasets/G4SAIDDATA2.0.tar.gz
wget https://geant4-data.web.cern.ch/datasets/G4ABLA3.1.tar.gz
wget https://geant4-data.web.cern.ch/datasets/G4INCL1.0.tar.gz
wget https://geant4-data.web.cern.ch/datasets/G4ENSDFSTATE2.3.tar.gz
tar -xzf *.tar.gz
ls -lah
cd $TOP


gsed -i 's@/Users/gcosmo/Software/release/install/@'$TOP'/Geant4-11.0.1-Darwin/@g'  $TOP/Geant4-11.0.1-Darwin/bin/geant4.sh

source $TOP/Geant4-11.0.1-Darwin/bin/geant4.sh
gsed -i 's@/opt/local/@/usr/local/Cellar/expat/2.4.7/@g'  Geant4-11.0.1-Darwin/lib/Geant4-11.0.1/Geant4PackageCache.cmake

mkdir -p test
cd test
cmake -B. -S ../examples/extended/eventgenerator  -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON -DCMAKE_INSTALL_PREFIX=$TOP/INSTALL -DGeant4_DIR=$TOP/Geant4-11.0.1-Darwin/lib/Geant4-11.0.1  -DPYTHIA6_ROOT_DIR=/usr/local -DHEPMC_DIR=/usr/local -DCMAKE_Fortran_FLAGS=-fPIC -DCMAKE_C_FLAGS=-fPIC -DCMAKE_CXX_FLAGS=-fPIC
cmake --build . -j 2
cmake --install .


#Running tests

cd  $TOP/test/HepMC/HepMCEx01  
$TOP/INSTALL/bin/HepMCEx01 hepmc_pygen.in
$TOP/INSTALL/bin/HepMCEx01 hepmc_ascii.in

cd  $TOP/test/HepMC/HepMCEx02
$TOP/INSTALL/bin/HepMCEx02 hepmc_pygen.in
$TOP/INSTALL/bin/HepMCEx02 hepmc_ascii.in
