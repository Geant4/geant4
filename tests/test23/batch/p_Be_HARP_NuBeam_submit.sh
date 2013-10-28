#!/usr/bin/env bash
#

APP_ENV=/home/yarba_j/work/test-geant4.9.6-ref08/geant4/tests/test23/batch/g4setup.sh

source ${APP_ENV}

# printenv | grep G4

PBS_WORK_DIR=${G4WORKDIR}/t23-bld

cd ${PBS_WORK_DIR}

echo "Executing Job " ${PBS_ARRAYID}

# uname -n

# JobID=1
JobID=${PBS_ARRAYID}
seed=$((1234+${JobID}))

printf "#verbose \n-1 \n#rad \n" >> NuBeam.Be.harp.${JobID}
printf "#events \n31250 \n" >> NuBeam.Be.harp.${JobID}
printf "#randomSeed\n" >> NuBeam.Be.harp.${JobID}
printf "%d\n" ${seed}  >> NuBeam.Be.harp.${JobID}
printf "#jobID\n" >> NuBeam.Be.harp.${JobID}
printf "%d\n" ${JobID} >> NuBeam.Be.harp.${JobID}
printf "//--------Proton_processes \n" >> NuBeam.Be.harp.${JobID}
printf "#particle \nproton \n#isHARP \n#position(mm) \n0. 0. -100. \n#direction \n0. 0. 1. \n" >> NuBeam.Be.harp.${JobID}

printf "//-------- \n" >> NuBeam.Be.harp.${JobID}

printf "#target-geom(mm) \n0.0 3.15 160. G4Tubs \n" >> NuBeam.Be.harp.${JobID}
printf "#material \nBe \n" >> NuBeam.Be.harp.${JobID}

printf "#momentum(MeV/c) \n8900.0 \n" >> NuBeam.Be.harp.${JobID}

printf "// --- \n#physicslist \nNuBeam \n#run \n" >> NuBeam.Be.harp.${JobID}

printf "#exit\n" >> NuBeam.Be.harp.${JobID}


./test23 NuBeam.Be.harp.${JobID}
#${PBS_WORK_DIR}/test23 NuBeam.Be.harp.${JobID}
# rm -f NuBeam.Be.harp.${JobID}
