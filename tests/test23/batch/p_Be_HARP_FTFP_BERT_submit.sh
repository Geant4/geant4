#!/usr/bin/env bash
#

APP_ENV=/home/yarba_j/work/test-geant4.9.6-ref08/geant4/tests/test23/batch/g4setup.sh

source ${APP_ENV}

PBS_WORK_DIR=${G4WORKDIR}/t23-bld

cd ${PBS_WORK_DIR}

# JobID=1
JobID=${PBS_ARRAYID}
seed=$((1234+${JobID}))

printf "#verbose \n-1 \n#rad \n" >> FTFP_BERT.Be.harp.${JobID}
printf "#events \n31250 \n" >> FTFP_BERT.Be.harp.${JobID}
printf "#randomSeed\n" >> FTFP_BERT.Be.harp.${JobID}
printf "%d\n" ${seed}  >> FTFP_BERT.Be.harp.${JobID}
printf "#jobID\n" >> FTFP_BERT.Be.harp.${JobID}
printf "%d\n" ${JobID} >> FTFP_BERT.Be.harp.${JobID}
printf "//--------Proton_processes \n" >> FTFP_BERT.Be.harp.${JobID}
printf "#particle \nproton \n#isHARP \n#position(mm) \n0. 0. -100. \n#direction \n0. 0. 1. \n" >> FTFP_BERT.Be.harp.${JobID}

printf "//-------- \n" >> FTFP_BERT.Be.harp.${JobID}

printf "#target-geom(mm) \n0.0 3.15 160. G4Tubs \n" >> FTFP_BERT.Be.harp.${JobID}
printf "#material \nBe \n" >> FTFP_BERT.Be.harp.${JobID}

printf "#momentum(MeV/c) \n8900.0 \n" >> FTFP_BERT.Be.harp.${JobID}

printf "// --- \n#physicslist \nftfp_bert \n#run \n" >> FTFP_BERT.Be.harp.${JobID}

printf "#exit\n" >> FTFP_BERT.Be.harp.${JobID}


./test23 FTFP_BERT.Be.harp.${JobID}
# ${PBS_WORK_DIR}/test23 FTFP_BERT.Be.harp.${JobID}
# rm -e -f FTFP_BERT.Be.harp.${JobID}
