#
# g4system.U
#
g4system="Linux-g++"

#
# g4dirs.U
#
g4install="/afs/cern.ch/user/s/sadilov/geant4"
g4include=""
g4base="/afs/cern.ch/user/s/sadilov/geant4/source"
g4workdir="/afs/cern.ch/user/s/sadilov/geant4"
g4tmp="/afs/cern.ch/user/s/sadilov/geant4/tmp"
g4lib="/afs/cern.ch/user/s/sadilov/geant4/lib"
g4bin="/afs/cern.ch/user/s/sadilov/geant4/bin"
g4data="/afs/cern.ch/user/s/sadilov/geant4/data"
g4levelgammadata="/afs/cern.ch/user/s/sadilov/geant4/data/PhotonEvaporation"

#
# g4clhep.U
#
g4clhep_base_dir="/usr/local"
g4clhep_include_dir="/usr/local/include"
g4clhep_lib_dir="/usr/local/lib"
g4clhep_lib="CLHEP"

#
# g4ospace
#
g4use_ospace=""
g4ospace_base_dir=""

#
# g4granular
#
g4global=""
g4lib_use_granular=""

#
# g4debug
#
g4debug=""

#
# g4shared
#
g4lib_build_shared=""
g4lib_build_archive="1"
g4lib_use_shared=""
g4lib_use_archive="1"

######################################
export G4SYSTEM="Linux-g++"

export G4INSTALL="/afs/cern.ch/user/s/sadilov/geant4"
export G4INCLUDE=""
export G4BASE="/afs/cern.ch/user/s/sadilov/geant4/source"
export G4WORKDIR="/afs/cern.ch/user/s/sadilov/geant4"
export G4TMP="/afs/cern.ch/user/s/sadilov/geant4/tmp"
export G4LIB="/afs/cern.ch/user/s/sadilov/geant4/lib"
export G4BIN="/afs/cern.ch/user/s/sadilov/geant4/bin"
export G4DATA="/afs/cern.ch/user/s/sadilov/geant4/data"
export G4LEVELGAMMADATA="/afs/cern.ch/user/s/sadilov/geant4/data/PhotonEvaporation"

export G4CLHEP_BASE_DIR="/usr/local"
export G4CLHEP_INCLUDE_DIR="/usr/local/include"
export G4CLHEP_LIB_DIR="/usr/local/lib"
export G4CLHEP_LIB="CLHEP"

export G4USE_OSPACE=""
export G4OSPACE_BASE_DIR=""

export G4LIB_USE_GRANULAR=""

export G4DEBUG=""

export G4LIB_BUILD_SHARED=""
export G4LIB_BUILD_ARCHIVE="1"
export G4LIB_USE_SHARED=""
export G4LIB_USE_ARCHIVE="1"

#####################################################################
echo "On this machine the G4SYSTEM=Linux-g++"

echo "On this machine the G4INSTALL=/afs/cern.ch/user/s/sadilov/geant4"
echo "On this machine the G4INCLUDE="
echo "On this machine the G4BASE=/afs/cern.ch/user/s/sadilov/geant4/source"
echo "On this machine the G4WORKDIR=/afs/cern.ch/user/s/sadilov/geant4"
echo "On this machine the G4TMP=/afs/cern.ch/user/s/sadilov/geant4/tmp"
echo "On this machine the G4LIB=/afs/cern.ch/user/s/sadilov/geant4/lib"
echo "On this machine the G4BIN=/afs/cern.ch/user/s/sadilov/geant4/bin"
echo "On this machine the G4DATA=/afs/cern.ch/user/s/sadilov/geant4/data"
echo "On this machine the G4LEVELGAMMADATA=/afs/cern.ch/user/s/sadilov/geant4/data/PhotonEvaporation"

echo "On this machine the G4CLHEP_BASE_DIR=/usr/local"
echo "On this machine the G4CLHEP_INCLUDE_DIR=/usr/local/include"
echo "On this machine the G4CLHEP_LIB_DIR=/usr/local/lib"
echo "On this machine the G4CLHEP_LIB=CLHEP"

echo "On this machine the G4USE_OSPACE="
echo "On this machine the G4OSPACE_BASE_DIR="

echo "On this machine the G4LIB_USE_GRANULAR="

echo "On this machine the G4DEBUG="

echo "On this machine the G4LIB_BUILD_SHARED="
echo "On this machine the G4LIB_BUILD_ARCHIVE=1"
echo "On this machine the G4LIB_USE_SHARED="
echo "On this machine the G4LIB_USE_ARCHIVE=1"

###############################################################
echo ""
echo "OK, going to /afs/cern.ch/user/s/sadilov/geant4/source and start 'make'..."
echo ""
cd /afs/cern.ch/user/s/sadilov/geant4/source
echo /home/sadilov/archive/g4conf.cvs/g4conf/save
if [ X = Xglobal ] ; then
make global
else
make
fi
if [ X = X1 ] ; then
make includes
fi
