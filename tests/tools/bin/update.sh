#! /bin/sh
# $Id: update.sh,v 1.9 2000-05-03 13:14:37 stesting Exp $
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
  echo DELETING G4RunManager.cc SO IT PICKS UP NEW TAG NAME...
  rm -f $G4INSTALL/source/run/src/G4RunManager.cc
  echo GETTING NEW DIRECTORIES...
else
  echo "DOING NOTHING (cvs -n)..."
  echo INSPECTING NEW DIRECTORIES...
fi

cd $G4INSTALL
echo "RUNNING updt.sh IN $G4INSTALL"
echo "REDIRECT Bonsai sdb FILE INTO INPUT, E.G., update.sh [-n] < bonsai.sdb"
echo "REDIRECT OUTPUT TO WHERE YOU LIKE, E.G., update.sh [-n] < bonsai.sdb >& update.log"
tests/tools/bin/updt.sh
