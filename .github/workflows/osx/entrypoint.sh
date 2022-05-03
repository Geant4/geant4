#!/bin/sh -l
set -x
export TOP=$(pwd)
mkdir LOCAL
cd LOCAL
brew tap davidchall/hep
brew install hepmc lhapdf wget

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

wget https://geant4-data.web.cern.ch/releases/lib_11.0.1/Darwin-clang13.0.0-Monterey.tar.gz

TOP=$(pwd)
mkdir -p test
cd test
cmake -B. -S ../examples/extended/eventgenerator  -DCMAKE_INSTALL_PREFIX=$TOP/INSTALL -DGeant4_DIR=/usr/lib64/Geant4-11.0.0  -DPYTHIA6_ROOT_DIR=/usr -DHEPMC_DIR=/usr -DCMAKE_Fortran_FLAGS=-fPIC -DCMAKE_C_FLAGS=-fPIC -DCMAKE_CXX_FLAGS=-fPIC
cmake --build . -j 2
cmake --install .
