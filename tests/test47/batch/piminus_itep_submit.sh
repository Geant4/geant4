#!/usr/bin/env bash
#

APP_ENV=/home/yarba_j/work/test-geant4.10.00-seq/geant4.10.00.b01/tests/test47/batch/g4setup.sh

source ${APP_ENV}

PBS_WORK_DIR=${G4WORKDIR}/test47

cd ${PBS_WORK_DIR}

# JobID=1
JobID=${PBS_ARRAYID}
seed=$((1234+${JobID}))

printf "#verbose \n0 \n#rad \n" >> ITEP.piminus.${JobID}
printf "#events \n31250 \n" >> ITEP.piminus.${JobID}
printf "#randomSeed\n" >> ITEP.piminus.${JobID}
printf "%d\n" ${seed}  >> ITEP.piminus.${JobID}
printf "#jobID\n" >> ITEP.piminus.${JobID}
printf "%d\n" ${JobID} >> ITEP.piminus.${JobID}
printf "//--------Piplus_processes \n" >> ITEP.piminus.${JobID}
printf "#particle \npi- \n#isITEP \n#position(mm) \n0. 0. 0. \n#direction \n0. 0. 1. \n//-------- \n" >> ITEP.piminus.${JobID}

printf "//-------- \n" >> ITEP.piminus.${JobID}
printf "#material \nC \n" >> ITEP.piminus.${JobID}
printf "#momentum(MeV/c) \n1400.0 \n" >> ITEP.piminus.${JobID}
printf "// --- \n#generator \nbertini \n#run \n// --- \n#generator \nbinary \n#run \n" >> ITEP.piminus.${JobID}

printf "//-------- \n" >> ITEP.piminus.${JobID}
printf "#momentum(MeV/c) \n5000.0 \n" >> ITEP.piminus.${JobID}
printf "// --- \n#generator \nbertini \n#run \n// --- \n#generator \nftfp \n#run \n" >> ITEP.piminus.${JobID}

printf "//-------- \n" >> ITEP.piminus.${JobID}
printf "#material \nU \n" >> ITEP.piminus.${JobID}
printf "#momentum(MeV/c) \n1400.0 \n" >> ITEP.piminus.${JobID}
printf "// --- \n#generator \nbertini \n#run \n// --- \n#generator \nbinary \n#run \n" >> ITEP.piminus.${JobID}

printf "//-------- \n" >> ITEP.piminus.${JobID}
printf "#momentum(MeV/c) \n5000.0 \n" >> ITEP.piminus.${JobID}
printf "// --- \n#generator \nbertini \n#run \n// --- \n#generator \nftfp \n#run \n" >> ITEP.piminus.${JobID}

printf "#exit\n" >> ITEP.piminus.${JobID}


${G4EXE}/test47 ITEP.piminus.${JobID}
