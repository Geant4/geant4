#!/bin/sh
#
#----------------------------------------------------------------------------
# This Bash shell script has the following 7 parameters:
#
#   1) Geant4 reference; e.g. 4.6.2.ref01 
#   2) Geant4 reference; e.g. 4.6.2.ref03
#   3) Physics List; e.g. LHEP
#   4) Calorimeter type; e.g. FeSci
#   5) Particle type; e.g. pi+
#   6) Beam Energy; e.g. 20GeV
#   7) Number of Events; e.g. 5k
#
# This script invokes the Python  writeSetup.py  that writes the setup
# and the Geant4 command files necessary to run the simulation of the 
# two Geant4 references. After that, this script invokes the Python
# script  dirStat/driver.py  that does the Statistical tests.
#----------------------------------------------------------------------------
#
echo ' ========== START simuDriver.sh ========== '
#
export REF1=$1
export REF2=$2
export PHYSICS=$3
export CALORIMETER=$4
export PARTICLE=$5
export ENERGY=$6
export EVENTS=$7
#
echo ' REF1        =' $REF1
echo ' REF2        =' $REF2
echo ' PHYSICS     =' $PHYSICS
echo ' CALORIMETER =' $CALORIMETER
echo ' PARTICLE    =' $PARTICLE
echo ' ENERGY      =' $ENERGY
echo ' EVENTS      =' $EVENTS
#
#--- Run the first reference ---
#
( export REF=$REF1 ;
  export LABEL=$REF-$PHYSICS-$CALORIMETER-$PARTICLE-$ENERGY-$EVENTS ;
  python2 writeSetup.py $REF $PHYSICS $CALORIMETER $PARTICLE $ENERGY $EVENTS ;
  mv run.g4 run.g4-$LABEL
  mv setup.sh setup.sh-$LABEL
  . setup.sh-$LABEL ;
  echo '  '; echo 'G4INSTALL = ' $G4INSTALL; echo 'running REF = ' $REF ; echo '  '
  mainStatAccepTest-$REF-$PHYSICS run.g4-$LABEL > output.log-$LABEL 2>&1 ;
  mv ntuple.hbook ntuple.hbook-$LABEL )
#
#--- Run the second reference ---
#
( export REF=$REF2 ; 
  export LABEL=$REF-$PHYSICS-$CALORIMETER-$PARTICLE-$ENERGY-$EVENTS ;
  python2 writeSetup.py $REF $PHYSICS $CALORIMETER $PARTICLE $ENERGY $EVENTS ;
  mv run.g4 run.g4-$LABEL
  mv setup.sh setup.sh-$LABEL
  . setup.sh-$LABEL ;
  echo '  '; echo 'G4INSTALL = ' $G4INSTALL; echo 'running REF = ' $REF ; echo '  '
  mainStatAccepTest-$REF-$PHYSICS run.g4-$LABEL > output.log-$LABEL 2>&1 ;
  mv ntuple.hbook ntuple.hbook-$LABEL )
#
#--- Run the statistical tests ---
#
( export LABEL=$PHYSICS-$CALORIMETER-$PARTICLE-$ENERGY-$EVENTS
  . setup.sh-$REF1-$LABEL
  cd dirStat/
  python2.2 driver.py $REF1 $REF2 $LABEL )
#
echo ' ========== END simuDriver.sh ========== '
#