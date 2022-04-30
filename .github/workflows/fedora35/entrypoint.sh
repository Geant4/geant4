#!/bin/sh -l
set -x
uname -a 
cat /etc/issue
yum -y install dnf-plugins-core
yum -y install gcc gcc-c++ gcc-gfortran make which cmake cmake-data cmake-filesystem
yum -y install HepMC3 HepMC3-devel HepMC HepMC-devel redhat-rpm-config
yum -y install expat-devel  xerces-c-devel xerces-c expat expat-devel zlib zlib-devel lhapdf lhapdf-devel
yum -y install pythia8-devel pythia8 pythia8-data
yum -y install yum-plugin-copr
yum -y copr enable averbyts/HEPrpms
yum -y install geant4 geant4-devel  clhep clhep-devel PTL-devel  pythia6 
yum -y install clean all
#Needed for tests, large
yum -y install geant4-data
yum -y install clean all

export FC=gfortran
export FCFLAGS=
TOP=$(pwd)
mkdir -p test
cd test
cmake -B. -S ../examples/extended/eventgenerator  -DCMAKE_INSTALL_PREFIX=$TOP/INSTALL -DGeant4_DIR=/usr/lib64/Geant4-11.0.0  -DPYTHIA6_ROOT_DIR=/usr -DHEPMC_DIR=/usr -DCMAKE_Fortran_FLAGS=-fPIC -DCMAKE_C_FLAGS=-fPIC -DCMAKE_CXX_FLAGS=-fPIC
cmake --build . -j 2
cmake --install .

#Running tests
source /usr/bin/geant4.sh

cd  $TOP/test/HepMC/HepMCEx01  
$TOP/INSTALL/bin/HepMCEx01 hepmc_pygen.in
$TOP/INSTALL/bin/HepMCEx01 hepmc_ascii.in

cd  $TOP/test/HepMC/HepMCEx02
$TOP/INSTALL/bin/HepMCEx02 hepmc_pygen.in
$TOP/INSTALL/bin/HepMCEx02 hepmc_ascii.in


out=$?
echo ::set-output name=out::$out
