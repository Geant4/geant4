#!/usr/bin/env bash
#

APP_ENV=/home/yarba_j/work/test-geant4.9.6-ref08/geant4/tests/test23/batch/g4setup.sh

source ${APP_ENV}

PBS_WORK_DIR=${G4WORKDIR}/t23-bld

cd ${PBS_WORK_DIR}

# JobID=1
JobID=${PBS_ARRAYID}
seed=$((1234+${JobID}))

printf "#verbose \n-1 \n#rad \n" >> QGSP_BERT.Be.harp.${JobID}
printf "#events \n31250 \n" >> QGSP_BERT.Be.harp.${JobID}
printf "#randomSeed\n" >> QGSP_BERT.Be.harp.${JobID}
printf "%d\n" ${seed}  >> QGSP_BERT.Be.harp.${JobID}
printf "#jobID\n" >> QGSP_BERT.Be.harp.${JobID}
printf "%d\n" ${JobID} >> QGSP_BERT.Be.harp.${JobID}
printf "//--------Proton_processes \n" >> QGSP_BERT.Be.harp.${JobID}
printf "#particle \nproton \n#isHARP \n#position(mm) \n0. 0. -100. \n#direction \n0. 0. 1. \n" >> QGSP_BERT.Be.harp.${JobID}

printf "//-------- \n" >> QGSP_BERT.Be.harp.${JobID}

printf "#target-geom(mm) \n0.0 3.15 160. G4Tubs \n" >> QGSP_BERT.Be.harp.${JobID}
printf "#material \nBe \n" >> QGSP_BERT.Be.harp.${JobID}

printf "#momentum(MeV/c) \n8900.0 \n" >> QGSP_BERT.Be.harp.${JobID}

printf "// --- \n#physicslist \nqgsp_bert \n#run \n" >> QGSP_BERT.Be.harp.${JobID}

printf "#exit\n" >> QGSP_BERT.Be.harp.${JobID}


./test23 QGSP_BERT.Be.harp.${JobID}
# ${PBS_WORK_DIR}/test23 QGSP_BERT.Be.harp.${JobID}
# rm -e -f QGSP_BERT.Be.harp.${JobID}
