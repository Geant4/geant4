#!/bin/bash
#
# Usage:
# process_examples.sh script 

#set -x

CURDIR=`pwd`

SCRIPT=$1
SCRIPT_ARG=$2

# call a script 
# {1} example directory
call_script() { 
  echo "... calling ${1} for example ${2}"
  ${1}
}  

# extended
#CATEGORY errorpropagation; polarisation
for DIR in extended; do
  echo ... processing $DIR
  cd $DIR
  EXTENDED_DIR=`pwd`
  for CATEGORY in `ls`; do
    # select directories with .README 
    if [ -f ${CATEGORY}/.README ]; then 
      echo ... processing $CATEGORY
      cd $CATEGORY
      CATEGORY_DIR=`pwd`
      for FILE in `find -name .README`; do
        EXAMPLE_DIR=`echo $FILE | sed sY/.READMEYYg`
        NOT_EXAMPLE=`echo ". ./gdml ./pythia ./MPI" | grep $EXAMPLE_DIR`
        if [ "${NOT_EXAMPLE}" = "" ]; then
          cd $EXAMPLE_DIR
          call_script $SCRIPT $EXAMPLE_DIR
          cd $CATEGORY_DIR
        fi
      done
    fi
    cd $EXTENDED_DIR
  done  
  # categories with a single example in the same level
  for EXAMPLE_DIR in errorpropagation; do
    cd $EXAMPLE_DIR
    call_script $SCRIPT $EXAMPLE_DIR
    cd $CATEGORY_DIR
    cd $EXTENDED_DIR
  done  
  cd $CURDIR
done         
