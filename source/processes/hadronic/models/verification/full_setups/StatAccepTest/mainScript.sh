#!/bin/sh
#
#----------------------------------------------------------------------------
# This Bash-shell script is the main, starting one.
# In  ***LOOKHERE***  you have to specify which Geant4 references you want
# to use, which Physics List, which calorimeter setup, which beam particle,
# which beam energy, which number of events.
# This program assumes the existence of the following executables:
#     mainStatAccepTest-$REF-$PHYSICS
# for the Geant4 reference and Physics List considered.
# This script invokes the following Python scripts:
#     writeSetup.py
#----------------------------------------------------------------------------
#
echo '========== START mainScript.sh ========== '
#
#***LOOKHERE***
export REF1=4.6.2.ref01
export REF2=4.6.2.ref03
export PHYSICS=LHEP
export CALORIMETER=FeSci
export PARTICLE=p
export ENERGY=20GeV
export EVENTS=5k
#export EVENTS=100
#***endLOOKHERE***
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
export LABEL=$PHYSICS-$CALORIMETER-$PARTICLE-$ENERGY-$EVENTS
. setup.sh-$REF1-$LABEL
cd dirStat/
python2.2 driver.py $REF1 $REF2 $LABEL
#
echo '========== END mainScript.sh ========== '
#