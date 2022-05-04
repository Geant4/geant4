#!/bin/sh -l
set -x
export TOP=$(pwd)
brew tap davidchall/hep
brew install lhapdf wget cmake coreutils  
brew install gnu-sed
brew install gcc
brew install expat
brew install pythia
brew install libxmu

cat <<EOF > hepmc2.rb
class Hepmc2 < Formula
  desc "C++ event record for Monte Carlo generators"
  homepage "https://hepmc.web.cern.ch/"
  url "https://gitlab.cern.ch/hepmc/HepMC/-/archive/2.06.11/HepMC-2.06.11.tar.gz"
  sha256 "ceaced62d39e4e2a1469fa2f20662d4d370279b3209930250766db02f44ae8de"

  option "with-test", "Test during installation"

  depends_on "cmake" => :build

  def install
    mkdir "../build" do
      system "cmake", buildpath, "-Dmomentum:STRING=GEV", "-Dlength:STRING=MM", *std_cmake_args
      system "make"
      system "make", "test" if build.with? "test"
      system "make", "install"
    end
  end

  test do
    cp_r share/"HepMC/examples/.", testpath
    system "make", "example_BuildEventFromScratch.exe"
    system "./example_BuildEventFromScratch.exe"
  end
end
EOF

brew install --build-from-source ./hepmc2.rb

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

mkdir -p Geant4-11.0.1-Darwin/share/Geant4-11.0.1/data/
cd Geant4-11.0.1-Darwin/share/Geant4-11.0.1/data/
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
