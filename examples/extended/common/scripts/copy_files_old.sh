#!/bin/sh

# Script to copy classes from common repository in examples.
# Usage:
# copy_files.sh className versionNumber [messenger] 

COMMON_DIR=../../../common
PREFIX="ExG4"
TO=`pwd`

fileName=$1
fileVersion=$2
isMessenger=$3

fileDir=""
if [ "$fileName" = "HbookAnalysisManager" ]; then
  fileDir="analysis"
  fileVersion=""
fi  
if [ "$fileName" = "DetectorConstruction" ]; then
  fileDir="detectorConstruction"
fi  
if [ "$fileName" = "PrimaryGeneratorAction" ]; then
  fileDir="primaryGenerator"
fi  
if [ "$fileName" = "EventAction" -o "$fileName" = "RunAction" ]; then
  fileDir="userActions"
fi  

echo "... copying $PREFIX$fileName$fileVersion"
cp $COMMON_DIR/$fileDir/include/$PREFIX$fileName$fileVersion".hh" include
cp $COMMON_DIR/$fileDir/src/$PREFIX$fileName$fileVersion".cc" src

if [ "$isMessenger" != "" ]; then
  echo "... copying $PREFIX$fileName"Messenger"$fileVersion"
  cp $COMMON_DIR/$fileDir/include/$PREFIX$fileName$fileVersion"Messenger.hh" include
  cp $COMMON_DIR/$fileDir/src/$PREFIX$fileName$fileVersion"Messenger.cc" src
fi
 
 
