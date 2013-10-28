#!/usr/bin/env bash
#

APP_ENV=/home/yarba_j/work/test-geant4.9.6-ref08/geant4/tests/test23/batch/g4setup.sh

source ${APP_ENV}

PBS_WORK_DIR=${G4WORKDIR}/t23-bld

cd ${PBS_WORK_DIR}

# JobID=1
JobID=${PBS_ARRAYID}
seed=$((1234+${JobID}))

printf "#verbose \n-1 \n#rad \n" >> NuBeam.C.na49.${JobID}
printf "#events \n31250 \n" >> NuBeam.C.na49.${JobID}
printf "#randomSeed\n" >> NuBeam.C.na49.${JobID}
printf "%d\n" ${seed}  >> NuBeam.C.na49.${JobID}
printf "#jobID\n" >> NuBeam.C.na49.${JobID}
printf "%d\n" ${JobID} >> NuBeam.C.na49.${JobID}
printf "//--------Proton_processes \n" >> NuBeam.C.na49.${JobID}
printf "#particle \nproton \n#isNA49 \n#position(mm) \n0. 0. -100. \n#direction \n0. 0. 1. \n" >> NuBeam.C.na49.${JobID}

printf "//-------- \n" >> NuBeam.C.na49.${JobID}

printf "#target-geom(mm) \n0.0 3.15 160. G4Tubs \n" >> NuBeam.C.na49.${JobID}
printf "#material \nC \n" >> NuBeam.C.na49.${JobID}

printf "#momentum(MeV/c) \n158000.0 \n" >> NuBeam.C.na49.${JobID}

printf "// --- \n#physicslist \nNuBeam \n#run \n" >> NuBeam.C.na49.${JobID}

printf "#exit\n" >> NuBeam.C.na49.${JobID}


./test23 NuBeam.C.na49.${JobID}
#${PBS_WORK_DIR}/bin/test23 NuBeam.C.na49.${JobID}
# rm -e -f NuBeam.C.na49.${JobID}
