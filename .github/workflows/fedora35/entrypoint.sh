#!/bin/sh -l
set -x
uname -a 
cat /etc/issue
yum -y install dnf-plugins-core
yum -y install  gcc gcc-c++ gcc-gfortran make which cmake cmake-data cmake-filesystem
yum -y install  HepMC3 HepMC3-devel HepMC HepMC-devel redhat-rpm-config
yum -y install  expat-devel  xerces-c-devel xerces-c expat expat-devel zlib zlib-devel lhapdf lhapdf-devel
yum -y install yum-plugin-copr
yum -y copr enable averbyts/HEPrpms
yum -y install geant4 geant4-devel  clhep clhep-devel PTL-devel  pythia6 
#yum -y install pythia8-devel pythia8 pythia8-data 

export FC=gfortran
export FCFLAGS=
mkdir -p test
cd test
cmake -B. -S ../examples/extended/eventgenerator -DGeant4_DIR=/usr/lib64/Geant4-11.0.0  -DPythia6_DIR=/usr -DHEPMC_DIR=/usr -DCMAKE_Fortran_FLAGS=-fPIC -DCMAKE_C_FLAGS=-fPIC -DCMAKE_CXX_FLAGS=-fPIC
make -j 2


out=$?
echo ::set-output name=out::$out
