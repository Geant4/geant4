#!/usr/bin/env bash
#

APP_ENV=/home/yarba_j/work/test-geant4.10.00-seq/geant4.10.00.b01/tests/test47/batch/g4setup.sh

source ${APP_ENV}

PBS_WORK_DIR=${G4WORKDIR}/test47

cd ${PBS_WORK_DIR}

# JobID=1
JobID=${PBS_ARRAYID}
seed=$((1234+${JobID}))

printf "#verbose\n0\n#rad\n" >> ITEP.piplus.${JobID}
printf "#events \n31250 \n" >> ITEP.piplus.${JobID}
printf "#randomSeed\n" >> ITEP.piplus.${JobID}
printf "%d\n" ${seed}  >> ITEP.piplus.${JobID}
printf "#jobID\n" >> ITEP.piplus.${JobID}
printf "%d\n" ${JobID} >> ITEP.piplus.${JobID}
printf "//--------Piplus_processes \n" >> ITEP.piplus.${JobID}
printf "#particle \npi+ \n#isITEP \n#position(mm) \n0. 0. 0. \n#direction \n0. 0. 1. \n//-------- \n" >> ITEP.piplus.${JobID}

printf "//-------- \n" >> ITEP.piplus.${JobID}
printf "#material \nC \n" >> ITEP.piplus.${JobID}
printf "#momentum(MeV/c) \n1400.0 \n" >> ITEP.piplus.${JobID}
printf "// --- \n#generator \nbertini \n#run \n// --- \n#generator \nbinary \n#run \n" >> ITEP.piplus.${JobID}

printf "//-------- \n" >> ITEP.piplus.${JobID}
printf "#momentum(MeV/c) \n5000.0 \n" >> ITEP.piplus.${JobID}
printf "// --- \n#generator \nbertini \n#run \n// --- \n#generator \nftfp \n#run \n" >> ITEP.piplus.${JobID}

printf "//-------- \n" >> ITEP.piplus.${JobID}
printf "#material \nU \n" >> ITEP.piplus.${JobID}
printf "#momentum(MeV/c) \n1400.0 \n" >> ITEP.piplus.${JobID}
printf "// --- \n#generator \nbertini \n#run \n// --- \n#generator \nbinary \n#run \n" >> ITEP.piplus.${JobID}

printf "//-------- \n" >> ITEP.piplus.${JobID}
printf "#momentum(MeV/c) \n5000.0 \n" >> ITEP.piplus.${JobID}
printf "// --- \n#generator \nbertini \n#run \n// --- \n#generator \nftfp \n#run \n" >> ITEP.piplus.${JobID}

printf "#exit\n" >> ITEP.piplus.${JobID}


${G4EXE}/test47 ITEP.piplus.${JobID}
