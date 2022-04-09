#!/bin/sh -l
set -x
uname -a 
cat /etc/issue
yum -y install dnf-plugins-core
yum -y install  gcc gcc-c++ gcc-gfortran make which cmake cmake-data cmake-filesystem
yum -y install  HepMC3 HepMC3-devel HepMC HepMC-devel
yum -y install yum-plugin-copr
yum -y copr enable averbyts/HEPrpms
yum -y install geant4 geant4-devel  clhep clhep-devel PTL-devel expat expat-devel zlib zlib-devel pythia6 lhapdf lhapdf-devel
#yum -y install pythia8-devel pythia8 pythia8-data 

export FC=gfortran
export FCFLAGS=
mkdir -p test
cd test
cmake -B. -S ../examples/extended/eventgenerator -DGeant4_DIR=/usr/lib64/Geant4-10.7.0  -DHepMC3_DIR=/usr/share/HepMC3/cmake -DCMAKE_Fortran_FLAGS=-fPIC -DCMAKE_C_FLAGS=-fPIC -DCMAKE_CXX_FLAGS=-fPIC
make -j 2


out=$?
echo ::set-output name=out::$out
