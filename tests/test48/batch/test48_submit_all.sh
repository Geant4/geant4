#!/usr/bin/env bash
#

source /products/setup

if [ `echo ${PATH} | grep pbs` ]; then
echo "PBS is already set"
else
export PATH=${PATH}:/usr/local/pbs/bin
fi

# G4RELEASE=${1}

qsub -q grunt -v G4RELEASE=${1},G4REGRESSION=${2} test48_run_piminus.sh

# NOTE: as of Oct. 2014 (geant4-10-00-ref-08), validation of mu- capture is suspended
#       because usage of the new mu- capture model (G4MuonMinusCapture) may result
#       in a stuck job, sujcet to a choice of startup random seed; infinite loop suspected

# all targets and models in 1 job
#
#qsub -q grunt -v G4RELEASE=${1},G4REGRESSION=${2} test48_run_muminus_chain.sh

# separate targets; 2 models in each job
#
#qsub -q grunt -v G4RELEASE=${1},G4REGRESSION=${2},target=Al test48_run_muminus.sh
#qsub -q grunt -v G4RELEASE=${1},G4REGRESSION=${2},target=Si test48_run_muminus.sh
#qsub -q grunt -v G4RELEASE=${1},G4REGRESSION=${2},target=Ca test48_run_muminus.sh
#qsub -q grunt -v G4RELEASE=${1},G4REGRESSION=${2},target=Fe test48_run_muminus.sh
#qsub -q grunt -v G4RELEASE=${1},G4REGRESSION=${2},target=Ag test48_run_muminus.sh
#qsub -q grunt -v G4RELEASE=${1},G4REGRESSION=${2},target=I  test48_run_muminus.sh
#qsub -q grunt -v G4RELEASE=${1},G4REGRESSION=${2},target=Au test48_run_muminus.sh
#qsub -q grunt -v G4RELEASE=${1},G4REGRESSION=${2},target=Pb test48_run_muminus.sh
#qsub -q grunt -v G4RELEASE=${1},G4REGRESSION=${2},target=S  test48_run_muminus.sh

#test_list=`ls -l | grep run | awk '{print $9}'`
#for test in ${test_list} ; do
#  echo "... Processing ... ${test} ... for ${1}"
##   qsub -q grunt -v ${1} ${test}
#  qsub -l nodes=1,walltime=48:00:00 -q grunt -v G4RELEASE=${1},G4REGRESSION=${2} ${test}
#done
