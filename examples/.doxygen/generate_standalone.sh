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
  #echo "processing ${1}"
  cd ../extended/${1}
  EXAMPLE_NAME=${PWD##*/}
  #echo "Generating Doxyfile for example: ${EXAMPLE_NAME}"
  cat $CURDIR/Doxyfile_standalone | sed s/"EXAMPLE_NAME"/"${EXAMPLE_NAME}"/g | sed sY"BACK_PATH"Y"${2}"Yg | sed sY"ADD_SHARED"Y"${ADD_SHARED}"Yg | sed sY"ADD_COMMON"Y"${ADD_COMMON}"Yg > Doxyfile
  ln -s $CURDIR/geant4.tag geant4.tag
  doxygen >& $CURDIR/doxygen_${EXAMPLE_NAME}.out
  rm Doxyfile
  rm geant4.tag
  cd $CURDIR
}

# process examples in second level directory in extended
for DIR in analysis/AnaEx01 analysis/AnaEx02 analysis/AnaEx03 electromagnetic/TestEm0 electromagnetic/TestEm1 electromagnetic/TestEm2 electromagnetic/TestEm3 electromagnetic/TestEm4 electromagnetic/TestEm5 electromagnetic/TestEm6 electromagnetic/TestEm7 electromagnetic/TestEm8 electromagnetic/TestEm9 electromagnetic/TestEm10 electromagnetic/TestEm11 electromagnetic/TestEm12 electromagnetic/TestEm13 electromagnetic/TestEm14 electromagnetic/TestEm15 electromagnetic/TestEm16 electromagnetic/TestEm17 electromagnetic/TestEm18 eventgenerator/particleGun eventgenerator/exgps eventgenerator/userPrimaryGenerator exoticphysics/channeling  exoticphysics/dmparticle exoticphysics/monopole geometry/transforms hadronic/Hadr00 hadronic/Hadr01 hadronic/Hadr02 hadronic/Hadr03 hadronic/Hadr04 hadronic/Hadr05 hadronic/Hadr06 hadronic/Hadr07 hadronic/Hadr08 hadronic/Hadr09 hadronic/Hadr10 hadronic/NeutronSource medical/electronScattering medical/electronScattering2 medical/fanoCavity medical/fanoCavity2 medical/GammaTherapy optical/OpNovice2 physicslists/extensibleFactory physicslists/factory physicslists/genericPL polarisation/Pol01 radioactivedecay/Activation radioactivedecay/rdecay01 radioactivedecay/rdecay02 runAndEvent/RE07; do
  generate ${DIR} ${BACK_PATH2}
done  

# process examples in third level directory in extended
for DIR in eventgenerator/HepMC/MCTruth eventgenerator/pythia/py8decayer hadronic/ParticleFluence/Calo hadronic/ParticleFluence/ConcentricSpheres hadronic/ParticleFluence/Layer hadronic/ParticleFluence/Sphere medical/dna/chem1 medical/dna/chem2 medical/dna/chem3 medical/dna/chem4 medical/dna/chem5 medical/dna/chem6 medical/dna/clustering medical/dna/dnadamage1 medical/dna/dnaphysics medical/dna/icsd medical/dna/jetcounter medical/dna/mfp medical/dna/microdosimetry medical/dna/microprox medical/dna/microyz medical/dna/moleculardna medical/dna/neuron medical/dna/pdb4dna medical/dna/range medical/dna/slowing medical/dna/scavenger medical/dna/splitting medical/dna/spower medical/dna/svalue medical/dna/wholeNuclearDNA medical/dna/wvalue parameterisations/gflash/gflasha parallel/TopC/ParN02 parallel/TopC/ParN04; do
  generate ${DIR} ${BACK_PATH3}
done  

# process examples in fourth level directory in extended
for DIR in  parallel/MPI/examples/exMPI01 parallel/MPI/examples/exMPI02  parallel/MPI/examples/exMPI03; do
  generate ${DIR} ${BACK_PATH4}
done  
 
cd $CURDIR
