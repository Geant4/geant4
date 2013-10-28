#!/usr/bin/env bash
#

APP_ENV=/home/yarba_j/work/test-geant4.9.6-ref08/geant4/tests/test23/batch/g4setup.sh

source ${APP_ENV}

PBS_WORK_DIR=${G4WORKDIR}/t23-bld

cd ${PBS_WORK_DIR}

# JobID=1
JobID=${PBS_ARRAYID}
seed=$((1234+${JobID}))

printf "#verbose \n-1 \n#rad \n" >> FTFP_BERT.C.na49.${JobID}
printf "#events \n31250 \n" >> FTFP_BERT.C.na49.${JobID}
printf "#randomSeed\n" >> FTFP_BERT.C.na49.${JobID}
printf "%d\n" ${seed}  >> FTFP_BERT.C.na49.${JobID}
printf "#jobID\n" >> FTFP_BERT.C.na49.${JobID}
printf "%d\n" ${JobID} >> FTFP_BERT.C.na49.${JobID}
printf "//--------Proton_processes \n" >> FTFP_BERT.C.na49.${JobID}
printf "#particle \nproton \n#isNA49 \n#position(mm) \n0. 0. -100. \n#direction \n0. 0. 1. \n" >> FTFP_BERT.C.na49.${JobID}

printf "//-------- \n" >> FTFP_BERT.C.na49.${JobID}

printf "#target-geom(mm) \n0.0 3.15 160. G4Tubs \n" >> FTFP_BERT.C.na49.${JobID}
printf "#material \nC \n" >> FTFP_BERT.C.na49.${JobID}

printf "#momentum(MeV/c) \n158000.0 \n" >> FTFP_BERT.C.na49.${JobID}

printf "// --- \n#physicslist \nftfp_bert \n#run \n" >> FTFP_BERT.C.na49.${JobID}

printf "#exit\n" >> FTFP_BERT.C.na49.${JobID}


./test23 FTFP_BERT.C.na49.${JobID}
#${PBS_WORK_DIR}/bin/test23 FTFP_BERT.C.na49.${JobID}
# rm -e -f FTFP_BERT.C.na49.${JobID}
