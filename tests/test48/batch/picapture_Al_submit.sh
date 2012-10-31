#!/usr/bin/env bash
#

APP_ENV=/home-cluck/yarba_j/work/test-central-geant4.9.6.b01/tests/test48/batch/g4setup.sh

source ${APP_ENV}

PBS_WORK_DIR=${G4WORKDIR}

cd ${PBS_WORK_DIR}

# JobID=1
JobID=${PBS_ARRAYID}
seed=$((1234+${JobID}))

printf "#verbose \n0 \n#rad \n" >> piminus.Al.${JobID}
printf "#events \n31250 \n" >> piminus.Al.${JobID}
printf "#randomSeed\n" >> piminus.Al.${JobID}
printf "%d\n" ${seed}  >> piminus.Al.${JobID}
printf "#jobID\n" >> piminus.Al.${JobID}
printf "%d\n" ${JobID} >> piminus.Al.${JobID}
printf "//--------Piminus_processes \n" >> piminus.Al.${JobID}
printf "#particle \npi- \n#position(mm) \n0. 0. 0. \n#direction \n0. 0. 1. \n//-------- \n" >> piminus.Al.${JobID}

printf "#material \nAl \n" >> piminus.Al.${JobID}

printf "#momentum(MeV/c) \n0.0 \n" >> piminus.Al.${JobID}
printf "// --- \n#generator \nBertiniPreCo \n#run \n#generator \nCHIPS \n#run \n#generator \nstopping \n#run \n" >> piminus.Al.${JobID}

printf "#exit\n" >> piminus.Al.${JobID}


${G4EXE}/test48 piminus.Al.${JobID}
