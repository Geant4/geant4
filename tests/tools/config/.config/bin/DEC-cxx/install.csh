#
# g4system.U
#
set g4system="DEC-cxx"

#
# g4dirs.U
#
set g4install="/usr/users/sadilov/geant4"
set g4base="/usr/users/sadilov/geant4/source"
set g4workdir="/usr/users/sadilov/geant4"
set g4tmp="/usr/users/sadilov/geant4/tmp"
set g4lib="/usr/users/sadilov/geant4/lib"
set g4bin="/usr/users/sadilov/geant4/bin"
set g4data="/usr/users/sadilov/geant4/data"
set g4levelgammadata="/usr/users/sadilov/geant4/data/PhotonEvaporation"

#
# g4clhep.U
#
set g4clhep_base_dir="/usr/clhep"
set g4clhep_include_dir="/usr/local/include"
set g4clhep_lib_dir="/usr/local/lib"
set g4clhep_lib="CLHEP"

#
# g4ospace.U
#
set g4use_ospace="n"
set g4ospace_base_dir="/usr/users/sadilov/ObjectSpace"

#################################################
setenv G4SYSTEM "DEC-cxx"

setenv G4INSTALL "/usr/users/sadilov/geant4"
setenv G4BASE "/usr/users/sadilov/geant4/source"
setenv G4WORKDIR "/usr/users/sadilov/geant4"
setenv G4TMP "/usr/users/sadilov/geant4/tmp"
setenv G4LIB "/usr/users/sadilov/geant4/lib"
setenv G4BIN "/usr/users/sadilov/geant4/bin"
setenv G4DATA "/usr/users/sadilov/geant4/data"
setenv G4LEVELGAMMADATA "/usr/users/sadilov/geant4/data/PhotonEvaporation"

setenv G4CLHEP_BASE_DIR "/usr/clhep"
setenv G4CLHEP_INCLUDE_DIR "/usr/local/include"
setenv G4CLHEP_LIB_DIR "/usr/local/lib"
setenv G4CLHEP_LIB "CLHEP"

setenv G4USE_OSPACE "n"
setenv G4OSPACE_BASE_DIR "/usr/users/sadilov/ObjectSpace"

####################################################
echo "On this machine the G4SYSTEM=DEC-cxx"

echo "On this machine the G4INSTALL=/usr/users/sadilov/geant4"
echo "On this machine the G4BASE=/usr/users/sadilov/geant4/source"
echo "On this machine the G4WORKDIR=/usr/users/sadilov/geant4"
echo "On this machine the G4TMP=/usr/users/sadilov/geant4/tmp"
echo "On this machine the G4LIB=/usr/users/sadilov/geant4/lib"
echo "On this machine the G4BIN=/usr/users/sadilov/geant4/bin"
echo "On this machine the G4DATA=/usr/users/sadilov/geant4/data"
echo "On this machine the G4LEVELGAMMADATA=/usr/users/sadilov/geant4/data/PhotonEvaporation"

echo "On this machine the G4CLHEP_BASE_DIR=/usr/clhep"
echo "On this machine the G4CLHEP_INCLUDE_DIR=/usr/local/include"
echo "On this machine the G4CLHEP_LIB_DIR=/usr/local/lib"
echo "On this machine the G4CLHEP_LIB=CLHEP"

echo "On this machine the G4USE_OSPACE=n"
echo "On this machine the G4OSPACE_BASE_DIR=/usr/users/sadilov/ObjectSpace"
###############################################################
echo ""
echo "OK, going to /usr/users/sadilov/geant4/source and start 'make'..."
echo ""
cd /usr/users/sadilov/geant4/source
echo /usr/users/sadilov/work/g4conf/config/scripts
if [ Xn = Xglobal ] ; then
make global
else
make
fi
