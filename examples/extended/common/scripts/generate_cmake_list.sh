#!/bin/sh

# Script to generate CMakeList.tx file from templates.
# The following fields in the template files are replaced with 
# specific example entries:
# EXAMPLE_NAME
# EXAMPLE_PROGRAM_NAME
# EXAMPLE_SCRIPT_LIST
# EXTERN_PACKAGE 
#
# Usage:     
# generate_cmake.sh categoryName
#     categoryName = all - process all categories
#
# By I. Hrivnacova, IPN Orsay

CURDIR=`pwd`
#FIX_XT="-lXt"
FIX_XT=""

# generate documentation 
# {1} example directory
# {2} external package name
generate() { 
  echo "  processing ${1}"
  cd ../../${1}
  EXAMPLE_NAME=${PWD##*/}
  
  # Extract program name
  # Get only first .cc and print a warning if more programs exist
  EXAMPLE_PROGRAM_FILES=`ls *.cc`
  FILE_FOUND="no"
  MESSAGE="no"
  for FILE in $EXAMPLE_PROGRAM_FILES; do
    if [ "${FILE_FOUND}" = "yes" ]; then
      if [ "${MESSAGE}" = "no" ]; then
        echo "  ... more programs found !!!"
        MESSAGE="yes" 
      fi  
      # update program name only if equal EXAMPLE_NAME
      # otherwise keep the first program found
      TMP_EXAMPLE_PROGRAM_NAME=`echo $FILE | sed  's/\(.*\)\..*/\1/'`
      if [ "${TMP_EXAMPLE_PROGRAM_NAME}" = "${EXAMPLE_NAME}" ]; then
        EXAMPLE_PROGRAM_NAME=$TMP_EXAMPLE_PROGRAM_NAME
        break
      fi  
    else
      EXAMPLE_PROGRAM_NAME=`echo $FILE | sed  's/\(.*\)\..*/\1/'`
      FILE_FOUND="yes"
    fi  
  done
  
  # Extract use of visualization
  GREP_VIS=`cat "$EXAMPLE_PROGRAM_NAME"".cc" | grep G4VisExecutive`
  if [ "${GREP_VIS}" = "" ]; then
    NO_VIS="_no_vis"
  else
    NO_VIS=""
  fi   
  
  # Extract list of scripts:
  # all files with extension *.mac *.in, *.out
  EXAMPLE_SCRIPT_LIST0=`ls *.mac *.in *.out *.gdml 2> /dev/null`
  EXAMPLE_SCRIPT_LIST=""
  for SCRIPT in $EXAMPLE_SCRIPT_LIST0; do
    EXAMPLE_SCRIPT_LIST="$EXAMPLE_SCRIPT_LIST ""$SCRIPT"
  done
    
  # Select a template if extern package
  if [ "${2}" = "" ]; then
    WITH_EXTERN=""
  else
    WITH_EXTERN="_with_extern"
    EXTERN_PACKAGE=${2}
    EXTERN_PACKAGE_TO_UPPER=`echo $EXTERN_PACKAGE | sed 'y/abcdefghijklmnopqrstuvwxyz/ABCDEFGHIJKLMNOPQRSTUVWXYZ/'`
  fi   
  if [ "${3}" = "optional" ]; then
    OPTIONAL="_opt"
  else
    OPTIONAL=""
  fi   

  # Generate file from templates
  # Project
  cat $CURDIR/CMakeLists_template_1.txt | sed s/"EXAMPLE_NAME"/"${EXAMPLE_NAME}"/g > CMakeLists.txt
  # Use Geant4 
  cat $CURDIR/CMakeLists_template_2"${NO_VIS}".txt >> CMakeLists.txt
  # Find extern package
  cat $CURDIR/CMakeLists_template_3"${WITH_EXTERN}""${OPTIONAL}".txt | sed s/"EXAMPLE_NAME"/"${EXAMPLE_NAME}"/g | sed s/"EXTERN_PACKAGE_TO_UPPER"/"${EXTERN_PACKAGE_TO_UPPER}"/g | sed s/"EXTERN_PACKAGE"/"${EXTERN_PACKAGE}"/g >> CMakeLists.txt
  # Locate sources and headers
  cat $CURDIR/CMakeLists_template_4"${WITH_EXTERN}".txt | sed s/"EXTERN_PACKAGE_TO_UPPER"/"${EXTERN_PACKAGE_TO_UPPER}"/g >> CMakeLists.txt
  # Add executable and link libraries
  cat $CURDIR/CMakeLists_template_5"${WITH_EXTERN}".txt | sed s/"EXAMPLE_PROGRAM_NAME"/"${EXAMPLE_PROGRAM_NAME}"/g | sed s/"EXTERN_PACKAGE_TO_UPPER"/"${EXTERN_PACKAGE_TO_UPPER}"/g | sed s/"FIX_XT"/"${FIX_XT}"/g >> CMakeLists.txt
  # Copy scripts
  cat $CURDIR/CMakeLists_template_6.txt | sed s/"EXAMPLE_NAME"/"${EXAMPLE_NAME}"/g | sed s/"EXAMPLE_SCRIPT_LIST"/"${EXAMPLE_SCRIPT_LIST}"/g >> CMakeLists.txt
  # Add custom target 
  if [ ! "${EXAMPLE_NAME}" = "${EXAMPLE_PROGRAM_NAME}" ]; then
    cat $CURDIR/CMakeLists_template_7.txt | sed s/"EXAMPLE_NAME"/"${EXAMPLE_NAME}"/g | sed s/"EXAMPLE_PROGRAM_NAME"/"${EXAMPLE_PROGRAM_NAME}"/g >> CMakeLists.txt
  fi
  # Install executable 
    cat $CURDIR/CMakeLists_template_8.txt | sed s/"EXAMPLE_PROGRAM_NAME"/"${EXAMPLE_PROGRAM_NAME}"/g >> CMakeLists.txt
    
  cd $CURDIR
}

# Functions per examples/extened categories
#

generate_analysis() {
  generate analysis/A01 AIDA optional 
  generate analysis/AnaEx01 
  generate analysis/AnaEx02 ROOT
  generate analysis/AnaEx03 AIDA
  generate analysis/N03Con
}

generate_biasing() {
  generate biasing/B01
  generate biasing/B02 AIDA
  generate biasing/ReverseMC01
}  

generate_common() {
  generate common/analysis HBOOK
  generate common/detectorConstruction  
  generate common/physicsList
  generate common/primaryGenerator
  generate common/userActions
}

generate_electromagnetic() {
  for DIR in electromagnetic/TestEm0 electromagnetic/TestEm1 electromagnetic/TestEm2 electromagnetic/TestEm3 electromagnetic/TestEm4 electromagnetic/TestEm5 electromagnetic/TestEm6 electromagnetic/TestEm7 electromagnetic/TestEm8 electromagnetic/TestEm9 electromagnetic/TestEm10 electromagnetic/TestEm11 electromagnetic/TestEm12 electromagnetic/TestEm13 electromagnetic/TestEm14 electromagnetic/TestEm15 electromagnetic/TestEm16 electromagnetic/TestEm17 electromagnetic/TestEm18; do
  #for DIR in electromagnetic/TestEm1; do
    generate ${DIR}
  done  
}  

generate_errorpropagation() {
  generate errorpropagation
}  

generate_eventgenerator() {
  generate eventgenerator/exgps AIDA optional
  generate eventgenerator/HepMC/HepMCEx01 HepMC
  generate eventgenerator/HepMC/HepMCEx02 HepMC
  generate eventgenerator/HepMC/MCTruth HepMC
  generate eventgenerator/particleGun
  generate eventgenerator/pythia/decayer6 Pythia6 optional
}  

generate_exoticphysics() {
  generate exoticphysics/monopole
}  

generate_field() {
  generate field/BlineTracer
  generate field/field01
  generate field/field02
  generate field/field03
  generate field/field04
  generate field/field05
  generate field/field06
}  

generate_g3tog4() {
  generate g3tog4/clGeometry
}  

generate_geometry() {
  generate geometry/olap
  generate geometry/transforms
}  

generate_hadronic() {
  generate hadronic/Hadr00
  generate hadronic/Hadr01
  generate hadronic/Hadr02
  generate hadronic/Hadr03
}  

generate_medical() {
  generate medical/DICOM
  generate medical/electronScattering
  generate medical/electronScattering2
  generate medical/fanoCavity
  generate medical/fanoCavity2
  generate medical/GammaTherapy
}  

generate_optical() {
  generate optical/LXe
  generate optical/wls
}

generate_parameterisations() {
  generate parameterisations/gflash
}

generate_persistency() {
  generate persistency/P01 ROOT
  generate persistency/P02 ROOT
  generate persistency/P03
  generate persistency/gdml/G01
  generate persistency/gdml/G02
  generate persistency/gdml/G03
  generate persistency/gdml/G04
}  

generate_polarisation() {
  generate polarisation/Pol01 AIDA optional
}

generate_radioactivedecay() {
  generate radioactivedecay/rdecay01
  generate radioactivedecay/rdecay02
}

generate_runAndEvent() {
  generate runAndEvent/RE01
  generate runAndEvent/RE02
  generate runAndEvent/RE03
  generate runAndEvent/RE04
}  

generate_visualization() {
  generate visualization/perspective
  generate visualization/standalone
  generate visualization/userVisAction
}

# Check arguments
if [ $# -ne 1 ]; then
  echo " Usage:  "   
  echo " generate_cmake.sh categoryName"
  echo "     categoryName = all - process all categories "
  exit
fi  

# Process a selected category or all 
CATEGORY=$1
if [ "$CATEGORY" = "all" ]; then
  # to be done
  echo generate_analysis
  echo generate_electromagnetic
else
  echo "processing ${CATEGORY}"
  generate_${CATEGORY}
fi  
   
cd $CURDIR
