#!/bin/sh

# Script to generate standalone Doxygen documentation 
# for the examples with non unique class names.

CURDIR=`pwd`

BACK_PATH2="../../../.doxygen/doc"
BACK_PATH3="../../../../.doxygen/doc"
BACK_PATH4="../../../../../.doxygen/doc"

# generate documentation 
# {1} example directory
# {2} relative path from the example to the doc directory
generate() { 
  echo "processing ${1}"
  cd ../extended/${1}
  EXAMPLE_NAME=${PWD##*/}
  # add directory if shared
  if [ "${EXAMPLE_NAME}" = "shared" ]; then
    cd ..
    DIR_NAME=${PWD##*/}
    EXAMPLE_NAME=${DIR_NAME}"_shared"
    cd shared
  fi  
  # detect use of classes from shared
  ADD_SHARED=""
  if [ -f SharedFilesList.txt ]; then
    if [ "`grep  SHARED_CLASSES_LIST SharedFilesList.txt`" = "SHARED_CLASSES_LIST" ]; then
      ADD_SHARED="../shared"
    fi
  fi    
  # detect use of classes from common
  ADD_COMMON=""
  if [ -f SharedFilesList.txt ]; then
    if [ "`grep  COMMON_CLASSES_LIST SharedFilesList.txt`" = "COMMON_CLASSES_LIST" ]; then
      ADD_COMMON="../../common ../../../common"
    fi
  fi    
  #echo "Generating Doxyfile for example: ${EXAMPLE_NAME}"
  cat $CURDIR/Doxyfile_standalone | sed s/"EXAMPLE_NAME"/"${EXAMPLE_NAME}"/g | sed sY"BACK_PATH"Y"${2}"Yg | sed sY"ADD_SHARED"Y"${ADD_SHARED}"Yg | sed sY"ADD_COMMON"Y"${ADD_COMMON}"Yg > Doxyfile
  doxygen >& $CURDIR/doxygen_${EXAMPLE_NAME}.out
  rm Doxyfile
  cd $CURDIR
}

# process examples in second level directory in extended
for DIR in analysis/shared analysis/AnaEx01 analysis/AnaEx02 analysis/AnaEx03 electromagnetic/TestEm0 electromagnetic/TestEm1 electromagnetic/TestEm2 electromagnetic/TestEm3 electromagnetic/TestEm4 electromagnetic/TestEm5 electromagnetic/TestEm6 electromagnetic/TestEm7 electromagnetic/TestEm8 electromagnetic/TestEm9 electromagnetic/TestEm10 electromagnetic/TestEm11 electromagnetic/TestEm12 electromagnetic/TestEm13 electromagnetic/TestEm14 electromagnetic/TestEm15 electromagnetic/TestEm16 electromagnetic/TestEm17 electromagnetic/TestEm18 eventgenerator/particleGun eventgenerator/exgps eventgenerator/userPrimaryGenerator exoticphysics/channeling  exoticphysics/dmparticle exoticphysics/monopole geometry/transforms hadronic/Hadr00 hadronic/Hadr01 hadronic/Hadr02 hadronic/Hadr03 hadronic/Hadr04 hadronic/Hadr05 hadronic/Hadr06 hadronic/Hadr07 hadronic/NeutronSource medical/electronScattering medical/electronScattering2 medical/fanoCavity medical/fanoCavity2 medical/GammaTherapy polarisation/Pol01 radioactivedecay/Activation radioactivedecay/rdecay01 radioactivedecay/rdecay02; do
  generate ${DIR} ${BACK_PATH2}
done  

# process examples in third level directory in extended
for DIR in eventgenerator/HepMC/MCTruth  parallel/TopC/ParN02  parallel/TopC/ParN04 medical/dna/chem1 medical/dna/chem2 medical/dna/chem3 medical/dna/chem4 medical/dna/clustering medical/dna/dnaphysics medical/dna/microdosimetry medical/dna/microyz medical/dna/pdb4dna medical/dna/range medical/dna/svalue medical/dna/wholeNuclearDNA medical/dna/wvalue; do
  generate ${DIR} ${BACK_PATH3}
done  

# process examples in fourth level directory in extended
for DIR in  parallel/MPI/examples/exMPI01 parallel/MPI/examples/exMPI02  parallel/MPI/examples/exMPI03; do
  generate ${DIR} ${BACK_PATH4}
done  
 
cd $CURDIR
