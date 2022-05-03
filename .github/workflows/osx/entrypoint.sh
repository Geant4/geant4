#!/bin/sh -l
set -x
export TOP=$(pwd)
brew tap davidchall/hep
brew install hepmc lhapdf wget cmake coreutils  

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
tar -xzvf Darwin-clang13.0.0-Monterey.tar.gz
ls -lah

mkdir -p test
cd test
cmake -B. -S ../examples/extended/eventgenerator  -DCMAKE_INSTALL_PREFIX=$TOP/INSTALL -DGeant4_DIR=$TOP/Darwin-clang13.0.0-Monterey/Geant4-11.0.1-Darwin/lib64/Geant4-11.0.1  -DPYTHIA6_ROOT_DIR=/usr/local -DHEPMC_DIR=/usr/local -DCMAKE_Fortran_FLAGS=-fPIC -DCMAKE_C_FLAGS=-fPIC -DCMAKE_CXX_FLAGS=-fPIC
cmake --build . -j 2
cmake --install .
