#!/bin/sh

# Script to copy classes specified in SharedFilesList.txt file.
# Usage: copy_files.sh [root_dir example_source_dir]
#
# By I. Hrivnacova, IPN Orsay

IS_KEYWORD="no"
WHERE=""
EXAMPLE_PATH="undefined"
COMMON_CLASSES_PATH="undefined"
SHARED_CLASSES_PATH="undefined"

NEXT_EXAMPLE_PATH="no"
NEXT_COMMON_CLASSES_PATH="no"
NEXT_SHARED_CLASSES_PATH="no"
NEXT_COMMON_CLASSES_LIST="no"
NEXT_SHARED_CLASSES_LIST="no"

EXTH=".hh"
EXTC=".cc"

# Files directories
CURDIR=`pwd`
if [ "$#" = "0" ]; then  
  # default setting when script is called without arguments from the example directory
  # which is inside examples top directory and files are searched via relative paths
  ROOT_DIR=""
  EXAMPLE_SOURCE_DIR="."
  RELATIVE="yes"
else  
  if [ "$#" = "2" ]; then  
  # setting when script is called with arguments 
  # and files are searched via absolute paths
    ROOT_DIR=$1
    EXAMPLE_SOURCE_DIR=$2
    RELATIVE="no"
  else
    echo "Usage: copy_files.sh [root_dir example_source_dir]"
    exit 1
  fi  
fi

cd $EXAMPLE_SOURCE_DIR

# check if {1} is a keyword
is_keyword() {
  if [ "${1}" = "EXAMPLE" -o "${1}" = "COMMON_CLASSES_PATH" -o "${1}" = "SHARED_CLASSES_PATH" -o "${1}" = "COMMON_CLASSES_LIST" -o "${1}" = "SHARED_CLASSES_LIST" ]; then 
    IS_KEYWORD="yes"
  else  
    IS_KEYWORD="no"
  fi  
}  

# check if filename {1} is defined with an extension
# and if yes decide whether it should go in include (.hh and .icc)
# or is src (.cc) 
# pront error if extension is unkbown 
where() {
  EXTENSION="`echo ${1#*.}`"
  if [ $EXTENSION = "${1}" ]; then
    WHERE=""
  else 
    if [ $EXTENSION = "hh" -o  $EXTENSION = "icc" ]; then
      WHERE="include"
    else 
      if [ $EXTENSION = "cc" ]; then
        WHERE="src"
      else
        WHERE="nowhere"
      fi
    fi
  fi              
}  

# reset all NEXT_* variables
reset_all() {
  NEXT_EXAMPLE_PATH="no"
  NEXT_COMMON_CLASSES_PATH="no"
  NEXT_SHARED_CLASSES_PATH="no"
  NEXT_COMMON_CLASSES_LIST="no"
  NEXT_SHARED_CLASSES_LIST="no"
}

# copy the file {1} from {2} to {3} with keaword {4} in messages
copy_file() {
# parameters:
# {1} $WORD; {2} from_path; {3} to_path; {4} COMMON/SHARED
  echo "getting file ${1}"
  if [ ${2} = "undefined" ]; then
    echo "The path to ${3} classes (${4}_CLASSES_PATH) must be set first"
    exit 1
  fi 
  FILE=`find ${2} -name ${1}` 
  IS_FILE=`find . -name ${1}` 
  if [ "$FILE" = "" ]; then
    echo "${1} not found in ${2}"
  else
    if [ ! "$IS_FILE" = "" ]; then
      echo "${1} is already installed."
      echo "You have to run clean to install sources again"
      exit 1
    else
      cp $FILE ${3}
    fi  
  fi   
}

for WORD in `cat SharedFilesList.txt`; do 
  ### process keywords
  #echo $WORD
  is_keyword $WORD
  if [ "$IS_KEYWORD" = "yes" ]; then 
    ### reset keywords
    reset_all
    if [ $WORD = "EXAMPLE" ]; then
      NEXT_EXAMPLE_PATH="yes"
    fi  
    if [ $WORD = "COMMON_CLASSES_PATH" ]; then
      NEXT_COMMON_CLASSES_PATH="yes"
    fi  
    if [ $WORD = "SHARED_CLASSES_PATH" ]; then
      NEXT_SHARED_CLASSES_PATH="yes"
    fi  
    if [ $WORD = "COMMON_CLASSES_LIST" ]; then
      NEXT_COMMON_CLASSES_LIST="yes"
    fi  
    if [ $WORD = "SHARED_CLASSES_LIST" ]; then
      NEXT_SHARED_CLASSES_LIST="yes"
    fi  
  else
    # path to example
    if [ $NEXT_EXAMPLE_PATH = "yes" ]; then
      if [ $RELATIVE = "no" ]; then
        echo "Setting path to the example $WORD"
        EXAMPLE_PATH="/${WORD}/"
      else  
        EXAMPLE_PATH=""
      fi  
      NEXT_EXAMPLE_PATH="no"
    fi      
    # path to common
    if [ $NEXT_COMMON_CLASSES_PATH = "yes" ]; then
      echo "Setting path to common classes $WORD"
      COMMON_CLASSES_PATH=${ROOT_DIR}${EXAMPLE_PATH}${WORD}
      NEXT_COMMON_CLASSES_PATH="no"
    fi      
    # path to shared
    if [ $NEXT_SHARED_CLASSES_PATH = "yes" ]; then
      echo "Setting path to shared classes $WORD"
      SHARED_CLASSES_PATH=${ROOT_DIR}${EXAMPLE_PATH}${WORD}
      NEXT_SHARED_CLASSES_PATH="no"
    fi  
    # classes from common
    if [ $NEXT_COMMON_CLASSES_LIST = "yes" ]; then
      where $WORD
      if [ "$WHERE" != "" -a  "$WHERE" != "nowhere" ]; then
        copy_file $WORD $COMMON_CLASSES_PATH $WHERE "COMMON"
      else 
        copy_file $WORD$EXTH $COMMON_CLASSES_PATH include "COMMON"
        copy_file $WORD$EXTC $COMMON_CLASSES_PATH src "COMMON"
      fi  
    fi     
    # classes from shared
    if [ $NEXT_SHARED_CLASSES_LIST = "yes" ]; then
      where $WORD
      if [ "$WHERE" != "" -a  "$WHERE" != "nowhere" ]; then
        copy_file $WORD $SHARED_CLASSES_PATH $WHERE "SHARED"
      else 
        copy_file $WORD$EXTH $SHARED_CLASSES_PATH include "SHARED"
        copy_file $WORD$EXTC $SHARED_CLASSES_PATH src "SHARED"
      fi  
    fi 
  fi  
done  

cd $CURDIR

echo "end copy_files.sh"

