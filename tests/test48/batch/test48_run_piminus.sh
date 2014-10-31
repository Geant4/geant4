#!/usr/bin/env bash
#

source /products/setup
setup root v5_34_01 -q "e2:prof"

# setup G4 datasets
#
source /home/g4p/pbs/g4-had-validation/env-setup/g4-datasets-setup.sh

# setup workdir
#
if [ "x" == "x$G4WORKDIR" ] ; then
G4WORKDIR=${PBS_O_WORKDIR}/.. 
else
    echo "Variable says: $G4WORKDIR"
    echo "Variable PBS_O_WORKDIR says: $PBS_O_WORKDIR"
fi

cd ${G4WORKDIR}

# JobID=1
#JobID=${PBS_ARRAYID}
#seed=$((1234+${JobID}))

./test48 test48.piminus
./test48 test48.kminus
./test48 test48.sigmaminus
./test48 test48.antiproton

# $ROOTSYS/bin/root -b -p -q Plot.C\($index\)  
$ROOTSYS/bin/root -b -p -q PiMinusModels.C
$ROOTSYS/bin/root -b -p -q KMinusModels.C
$ROOTSYS/bin/root -b -p -q AntiProtonModels.C

if [ "x" == "x$G4REGRESSION" ]; then
   echo "Regression testing is NOT requested"
else
   echo "G4REGRESSION says: $G4REGRESSION"
$ROOTSYS/bin/root -b -p -q PiMinusRegre.C
$ROOTSYS/bin/root -b -p -q KMinusRegre.C
$ROOTSYS/bin/root -b -p -q AntiProtonRegre.C
fi
 
 
if [ "x" == "x$G4RELEASE" ] ; then
    echo "Variable G4RELEASE is not set"
else
#
if [ ! -d /home/g4p/pbs/g4-had-validation/results/${G4RELEASE} ]; then
mkdir /home/g4p/pbs/g4-had-validation/results/${G4RELEASE}
fi
if [ ! -d /home/g4p/pbs/g4-had-validation/results/${G4RELEASE}/test48 ]; then
mkdir /home/g4p/pbs/g4-had-validation/results/${G4RELEASE}/test48
fi
if [ ! -d /home/g4p/pbs/g4-had-validation/results/${G4RELEASE}/test48/models ]; then
mkdir /home/g4p/pbs/g4-had-validation/results/${G4RELEASE}/test48/models
fi
cp *models.gif /home/g4p/pbs/g4-had-validation/results/${G4RELEASE}/test48/models/.
#
if [ "x" == "x$G4REGRESSION" ]; then
   echo "NO regression plots"
else
# move regression plots to results area
if [ ! -d /home/g4p/pbs/g4-had-validation/results/${G4RELEASE}/test48/regre ]; then
mkdir /home/g4p/pbs/g4-had-validation/results/${G4RELEASE}/test48/regre
fi
cp *-regre.gif /home/g4p/pbs/g4-had-validation/results/${G4RELEASE}/test48/regre/.
fi
#
cp *.root /scratch-shared/g4p/pbs/g4-had-validation/regression-test-files/test48/${G4RELEASE}/.
#
fi
