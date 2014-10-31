#!/usr/bin/env bash
#

source /products/setup
setup root v5_34_01 -q "e2:prof"

# setup G4 datasets
#
source /home/g4p/pbs/g4-had-validation/env-setup/g4-datasets-setup.sh

# setup workdir
#
if [ "x" == "x$G4WORKDIR" ] ; then
G4WORKDIR=${PBS_O_WORKDIR}/.. 
else
    echo "Variable says: $G4WORKDIR"
    echo "Variable PBS_O_WORKDIR says: $PBS_O_WORKDIR"
fi

cd ${G4WORKDIR}

# JobID=1
JobID=${PBS_ARRAYID}
seed=$((1234+${JobID}))

printf "#verbose \n0 \n#rad \n" >> piminus.Pb.${JobID}
printf "#events \n31250 \n" >> piminus.Pb.${JobID}
printf "#randomSeed\n" >> piminus.Pb.${JobID}
printf "%d\n" ${seed}  >> piminus.Pb.${JobID}
printf "#jobID\n" >> piminus.Pb.${JobID}
printf "%d\n" ${JobID} >> piminus.Pb.${JobID}
printf "//--------Piminus_processes \n" >> piminus.Pb.${JobID}
printf "#particle \npi- \n#position(mm) \n0. 0. 0. \n#direction \n0. 0. 1. \n//-------- \n" >> piminus.Pb.${JobID}

printf "#material \nPb \n" >> piminus.Pb.${JobID}

printf "#momentum(MeV/c) \n0.0 \n" >> piminus.Pb.${JobID}
printf "// --- \n#generator \nBertiniPreCo \n#run \n#generator \nCHIPS \n#run \n#generator \nstopping \n#run \n" >> piminus.Pb.${JobID}

printf "#exit\n" >> piminus.Pb.${JobID}


${G4EXE}/test48 piminus.Pb.${JobID}
