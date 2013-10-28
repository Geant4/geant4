#!/usr/bin/env bash
#

APP_ENV=/home/yarba_j/work/test-geant4.9.6-ref08/geant4/tests/test23/batch/g4setup.sh

source ${APP_ENV}

PBS_WORK_DIR=${G4WORKDIR}/t23-bld

cd ${PBS_WORK_DIR}

# JobID=1
JobID=${PBS_ARRAYID}
seed=$((1234+${JobID}))

printf "#verbose \n-1 \n#rad \n" >> QGSP_BERT.C.na49.${JobID}
printf "#events \n31250 \n" >> QGSP_BERT.C.na49.${JobID}
printf "#randomSeed\n" >> QGSP_BERT.C.na49.${JobID}
printf "%d\n" ${seed}  >> QGSP_BERT.C.na49.${JobID}
printf "#jobID\n" >> QGSP_BERT.C.na49.${JobID}
printf "%d\n" ${JobID} >> QGSP_BERT.C.na49.${JobID}
printf "//--------Proton_processes \n" >> QGSP_BERT.C.na49.${JobID}
printf "#particle \nproton \n#isNA49 \n#position(mm) \n0. 0. -100. \n#direction \n0. 0. 1. \n" >> QGSP_BERT.C.na49.${JobID}

printf "//-------- \n" >> QGSP_BERT.C.na49.${JobID}

printf "#target-geom(mm) \n0.0 3.15 160. G4Tubs \n" >> QGSP_BERT.C.na49.${JobID}
printf "#material \nC \n" >> QGSP_BERT.C.na49.${JobID}

printf "#momentum(MeV/c) \n158000.0 \n" >> QGSP_BERT.C.na49.${JobID}

printf "// --- \n#physicslist \nqgsp_bert \n#run \n" >> QGSP_BERT.C.na49.${JobID}

printf "#exit\n" >> QGSP_BERT.C.na49.${JobID}


./test23 QGSP_BERT.C.na49.${JobID}
#${PBS_WORK_DIR}/bin/test23 QGSP_BERT.C.na49.${JobID}
# rm -e -f QGSP_BERT.C.na49.${JobID}
