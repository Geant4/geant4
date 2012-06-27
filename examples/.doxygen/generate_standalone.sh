#!/bin/sh

# Script to generate standalone Doxygen documentation 
# for the examples with non unique class names.

CURDIR=`pwd`

BACK_PATH2="../../../.doxygen/doc"
BACK_PATH3="../../../../.doxygen/doc"

# process examples in second level directory in extended
for DIR in analysis/shared analysis/AnaEx01 analysis/AnaEx02 analysis/AnaEx03 electromagnetic/TestEm0 electromagnetic/TestEm1 electromagnetic/TestEm2 electromagnetic/TestEm3 electromagnetic/TestEm4 electromagnetic/TestEm5 electromagnetic/TestEm6 electromagnetic/TestEm7 electromagnetic/TestEm8 electromagnetic/TestEm9 electromagnetic/TestEm10 electromagnetic/TestEm11 electromagnetic/TestEm12 electromagnetic/TestEm13 electromagnetic/TestEm14 electromagnetic/TestEm15 electromagnetic/TestEm16 electromagnetic/TestEm17 electromagnetic/TestEm18 eventgenerator/particleGun exoticphysics/monopole geometry/transforms hadronic/Hadr00 hadronic/Hadr01 hadronic/Hadr02 hadronic/Hadr03 medical/electronScattering medical/electronScattering2 medical/fanoCavity medical/fanoCavity2 medical/GammaTherapy parallel/ParN02 parallel/ParN04 polarisation/Pol01 radioactivedecay/rdecay01 radioactivedecay/rdecay02; do
#for DIR in  analysis/shared ; do
  echo "processing $DIR"
  cd ../extended/$DIR
  EXAMPLE_NAME=${PWD##*/}
  # add directory if shared
  if [ "${EXAMPLE_NAME}" = "shared" ]; then
    cd ..
    DIR_NAME=${PWD##*/}
    EXAMPLE_NAME=${DIR_NAME}"_shared"
    cd shared
  fi  
  #echo "Generating Doxyfile for example: ${EXAMPLE_NAME}"
  cat $CURDIR/Doxyfile_standalone | sed s/"EXAMPLE_NAME"/"${EXAMPLE_NAME}"/g | sed sY"BACK_PATH"Y"${BACK_PATH2}"Yg > Doxyfile
  doxygen >& $CURDIR/doxygen_${EXAMPLE_NAME}.out
  rm Doxyfile
  cd $CURDIR
done  

# process examples in third level directory in extended
for DIR in  parallel/MPI/exMPI01 parallel/MPI/exMPI02 persistency/gdml/G01 persistency/gdml/G02 persistency/gdml/G03; do
#for DIR in  parallel/MPI/exMPI02 ; do
  echo "processing $DIR"
  cd ../extended/$DIR
  EXAMPLE_NAME=${PWD##*/}
  #echo "Generating Doxyfile for example: ${EXAMPLE_NAME}"
  cat $CURDIR/Doxyfile_standalone | sed s/"EXAMPLE_NAME"/"${EXAMPLE_NAME}"/g | sed sY"BACK_PATH"Y"${BACK_PATH3}"Yg > Doxyfile
  doxygen >& $CURDIR/doxygen_${EXAMPLE_NAME}.out
  rm Doxyfile
  cd $CURDIR
done  
 
