#!/bin/sh
#
#----------------------------------------------------------------------------
# This Bash shell script has the following 10 parameters:
#
#   1)  Geant4 reference; e.g. 4.6.2.ref03 , or a local directory. 
#   2)  Flag to know whether the first reference should be simulated or not.
#   3)  Geant4 reference; e.g. 4.6.2.ref04 , or empty, or a local directory.
#   4)  Flag to know whether the second reference should be simulated or not.
#   5)  Flag to know whether the statistical tests should be run or not.
#   6)  Physics List; e.g. LHEP
#   7)  Calorimeter type; e.g. FeSci
#   8)  Particle type; e.g. pi+
#   9)  Beam Energy; e.g. 20GeV
#   10) Number of Events; e.g. 5k
#
# This script invokes the Python  build.py  that writes the setup,
# builds the executable, and writes the Geant4 command file.
# After that, this script eventually (if the flag is on) builds the
# executable  pvalue  and then calls the Python script  dirStat/driver.py  
# that does the Statistical tests.
#----------------------------------------------------------------------------
#
echo ' ========== START simuDriver.sh ========== '
#
export REF1=$1
export SIM_REF1=$2
export REF2=$3
export SIM_REF2=$4
export RUN_STAT=$5
export PHYSICS=$6
export CALORIMETER=$7
export PARTICLE=$8
export ENERGY=$9
export EVENTS=${10}
#
echo ' REF1        =' $REF1
echo ' SIM_REF1    =' $SIM_REF1
echo ' REF2        =' $REF2
echo ' SIM_REF2    =' $SIM_REF2
echo ' RUN_STAT    =' $RUN_STAT
echo ' PHYSICS     =' $PHYSICS
echo ' CALORIMETER =' $CALORIMETER
echo ' PARTICLE    =' $PARTICLE
echo ' ENERGY      =' $ENERGY
echo ' EVENTS      =' $EVENTS
#
#--- Run the first reference ---
#
( if [ X$SIM_REF1 == XYes ] ; then
    ###echo " I AM HERE 1 "
    export REF=$REF1 ;
    export LABEL=$REF-$PHYSICS-$CALORIMETER-$PARTICLE-$ENERGY-$EVENTS ;
    python2 build.py $REF $PHYSICS $CALORIMETER $PARTICLE $ENERGY $EVENTS ;
    mv run.g4 run.g4-$LABEL
    mv setup.sh setup.sh-$LABEL
    . setup.sh-$LABEL ;
    echo '  '; echo ' G4INSTALL = ' $G4INSTALL; echo ' running REF = ' $REF ; echo '  '
    rm -rf tmp/
    gmake
    mv bin/$G4SYSTEM/mainStatAccepTest bin/$G4SYSTEM/mainStatAccepTest-$REF-$PHYSICS
    mainStatAccepTest-$REF-$PHYSICS run.g4-$LABEL > output.log-$LABEL 2>&1 ;
    mv ntuple.hbook ntuple.hbook-$LABEL
  fi )
#
#--- Run the second reference ---
#
( if [ X$SIM_REF2 == XYes ] ; then
    ###echo " I AM HERE 2 "
    export REF=$REF2 ; 
    export LABEL=$REF-$PHYSICS-$CALORIMETER-$PARTICLE-$ENERGY-$EVENTS ;
    python2 build.py $REF $PHYSICS $CALORIMETER $PARTICLE $ENERGY $EVENTS ;
    mv run.g4 run.g4-$LABEL
    mv setup.sh setup.sh-$LABEL
    . setup.sh-$LABEL ;
    echo '  '; echo ' G4INSTALL = ' $G4INSTALL; echo ' running REF = ' $REF ; echo '  '
    rm -rf tmp/
    gmake
    mv bin/$G4SYSTEM/mainStatAccepTest bin/$G4SYSTEM/mainStatAccepTest-$REF-$PHYSICS
    mainStatAccepTest-$REF-$PHYSICS run.g4-$LABEL > output.log-$LABEL 2>&1 ;
    mv ntuple.hbook ntuple.hbook-$LABEL 
  fi )
#
#--- Run the statistical tests ---
#
( if [ X$RUN_STAT == XYes ] ; then
    ###echo " I AM HERE 3 "
    export LABEL=$PHYSICS-$CALORIMETER-$PARTICLE-$ENERGY-$EVENTS
    . setup.sh-$REF1-$LABEL
    cd dirStat/
    rm pvalue.o pvalue
    gmake
    python2.2 driver.py $REF1 $REF2 $LABEL 
  fi )
#
echo ' ========== END simuDriver.sh ========== '
#