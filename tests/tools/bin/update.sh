#! /bin/sh
# $Id: update.sh,v 1.7 2000-01-17 09:45:43 stesting Exp $
# Edit stt-prod.sdb or stt-dev.sdb and execute.
# Usage: update.sh [-n]

# Some checks :
if [ -z "${G4INSTALL}" ] ; then
  echo "G4INSTALL not set. Execute setup file first !"
  exit
fi
if [ -z "${G4SYSTEM}" ] ; then
  echo "G4SYSTEM not set. Execute setup file first !"
  exit
fi
if [ -z "${G4WORKDIR}" ] ; then
  echo "G4WORKDIR not set. Execute setup file first !"
  exit
fi

# Make $G4WORKDIR/stt directory :
dir=$G4WORKDIR/stt
if [ -d $dir ] ; then
 echo "" > /dev/null
else 
 echo "$dir does not exist."
 echo "This is now an error.  Create a directory or a symbolic link"
 echo "  to a directory by hand or in a calling script."
fi
#
# Make $G4WORKDIR/stt/$G4SYSTEM directory :
dir=$G4WORKDIR/stt/$G4SYSTEM
if [ -d $dir ] ; then
 echo "" > /dev/null
else 
 mkdir $dir
 echo "$dir created."
fi
#

for i in $*
do
case $i in
  -n) NOTHING=-n; export NOTHING;;
esac
done

echo

if [ X$NOTHING = X ]
then
  echo ACTUALLY UPDATING...
  echo DELETING G4RunManager.cc SO IT PICKS UP NEW '$Name: not supported by cvs2svn $'...
  rm -f $G4INSTALL/source/run/src/G4RunManager.cc
  echo GETTING NEW DIRECTORIES...
else
  echo "DOING NOTHING (cvs -n)..."
  echo INSPECTING NEW DIRECTORIES...
fi

cd $G4INSTALL/..
geant4/tests/tools/bin/updt.sh  < geant4/tests/stt-$REF.sdb \
  > $dir/update.log 2>&1
