#!/usr/bin/env bash
#

APP_ENV=/home/yarba_j/work/test-geant4.10.00-seq/geant4.10.00.b01/tests/test47/batch/g4setup.sh

source ${APP_ENV}

PBS_WORK_DIR=${G4WORKDIR}/test47

cd ${PBS_WORK_DIR}

# JobID=1
JobID=${PBS_ARRAYID}
seed=$((1234+${JobID}))

printf "#verbose\n0\n#rad\n" >> ITEP.proton.${JobID}
printf "#events \n31250 \n" >> ITEP.proton.${JobID}
printf "#randomSeed\n" >> ITEP.proton.${JobID}
printf "%d\n" ${seed}  >> ITEP.proton.${JobID}
printf "#jobID\n" >> ITEP.proton.${JobID}
printf "%d\n" ${JobID} >> ITEP.proton.${JobID}
printf "//--------Proton_processes \n" >> ITEP.proton.${JobID}
printf "#particle\nproton\n#isITEP\n#position(mm)\n0. 0. 0.\n#direction\n0. 0. 1.\n//-------- \n" >> ITEP.proton.${JobID}

printf "//-------- \n" >> ITEP.proton.${JobID}
printf "#material \nC \n" >> ITEP.proton.${JobID}
printf "#momentum(MeV/c) \n1400.0 \n" >> ITEP.proton.${JobID}
printf "// --- \n#generator \nbertini \n#run \n// --- \n#generator \nbinary \n#run \n" >> ITEP.proton.${JobID}

printf "//-------- \n" >> ITEP.proton.${JobID}
printf "#momentum(MeV/c) \n7500.0 \n" >> ITEP.proton.${JobID}
printf "// --- \n#generator \nbertini \n#run \n// --- \n#generator \nftfp \n#run \n" >> ITEP.proton.${JobID}

printf "//-------- \n" >> ITEP.proton.${JobID}
printf "#material \nU \n" >> ITEP.proton.${JobID}
printf "#momentum(MeV/c) \n1400.0 \n" >> ITEP.proton.${JobID}
printf "// --- \n#generator \nbertini \n#run \n// --- \n#generator \nbinary \n#run \n" >> ITEP.proton.${JobID}

printf "//-------- \n" >> ITEP.proton.${JobID}
printf "#momentum(MeV/c) \n7500.0 \n" >> ITEP.proton.${JobID}
printf "// --- \n#generator \nbertini \n#run \n// --- \n#generator \nftfp \n#run \n" >> ITEP.proton.${JobID}

printf "#exit\n" >> ITEP.proton.${JobID}


${G4EXE}/test47 ITEP.proton.${JobID}
