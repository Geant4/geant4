#!/bin/sh

# Script to clean classes specified in SharedFilesList.txt file.
# Usage: clean_files.sh [example_source_dir]
#
# By I. Hrivnacova, IPN Orsay

CURDIR=`pwd`

IS_KEYWORD="no"
WHERE=""
NEXT_COMMON_CLASSES_LIST="no"
NEXT_SHARED_CLASSES_LIST="no"

EXTH=".hh"
EXTC=".cc"

# Files directories
CURDIR=`pwd`
# default setting when script is called from the example directory
# which is inside examples top directory a files are searcxh via relative paths
if [ "$#" = "0" ]; then  
  EXAMPLE_SOURCE_DIR="."
else  
  if [ "$#" = "1" ]; then  
    EXAMPLE_SOURCE_DIR=$1
  else
    echo "Usage: clean_files.sh [example_source_dir]"
    exit 1
  fi  
fi

cd $EXAMPLE_SOURCE_DIR

# check if {1} is a keyword
is_keyword() {
  if [ "${1}" = "COMMON_CLASSES_PATH" -o "${1}" = "SHARED_CLASSES_PATH" -o "${1}" = "COMMON_CLASSES_LIST" -o "${1}" = "SHARED_CLASSES_LIST" ]; then 
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
  NEXT_COMMON_CLASSES_LIST="no"
  NEXT_SHARED_CLASSES_LIST="no"
}

# remove the file {1} 
remove_file() {
  echo "removing ${1}"
  FILE=`find . -name ${1}` 
  if [ "$FILE" = "" ]; then
    echo "${1} not found."
  else
    rm -f $FILE
  fi  
}

for WORD in `cat SharedFilesList.txt`; do 
  ### process keywords
  #echo $WORD
  is_keyword $WORD
  if [ "$IS_KEYWORD" = "yes" ]; then 
    ### reset keywords
    reset_all
    if [ $WORD = "COMMON_CLASSES_LIST" ]; then
      NEXT_COMMON_CLASSES_LIST="yes"
    fi  
    if [ $WORD = "SHARED_CLASSES_LIST" ]; then
      NEXT_SHARED_CLASSES_LIST="yes"
    fi  
  else
    # classes from common
    if [ $NEXT_COMMON_CLASSES_LIST = "yes" ]; then
      where $WORD
      if [ "$WHERE" != "" -a  "$WHERE" != "nowhere" ]; then
        remove_file $WORD
      else  
        remove_file $WORD$EXTH 
        remove_file $WORD$EXTC
      fi  
    fi     
    # classes from shared
    if [ $NEXT_SHARED_CLASSES_LIST = "yes" ]; then
      where $WORD
      if [ "$WHERE" != "" -a  "$WHERE" != "nowhere" ]; then
        remove_file $WORD
      else  
        remove_file $WORD$EXTH
        remove_file $WORD$EXTC
      fi  
    fi  
  fi  
done  
