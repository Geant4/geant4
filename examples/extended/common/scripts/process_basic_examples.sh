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

# basic
for DIR in basic; do
  echo ... processing $DIR
  cd $DIR
  BASIC_DIR=`pwd`
  for EXAMPLE in B1 B2/B2a B3 B4/B4a B4/B4b B4/B4c B4/B4d; do
    cd $EXAMPLE
    call_script $SCRIPT $EXAMPLE
    cd $BASIC_DIR
  done
  cd $CURDIR
done   
