#! /bin/sh
# $Id: update.sh,v 1.13 2000-08-01 08:24:56 stesting Exp $
# For tagset NNN, extract and check  bonsai<NNN>.sdb then run update.sh.
# Usage: update.sh [-n] < bonsai<NNN>.sdb >& update<NNN>.log

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
if [ -z "${G4STTDIR}" ] ; then
  echo "G4STTDIR not set. Execute setup file first !"
  exit
fi

# Assuming standard practice of running in tests/tools/bin
# ensure the environment variables match the current directory
# /afs/cern.ch/sw/geant4/stt/dev1/src/geant4

mysrc=`echo $G4INSTALL | cut -f 7 -d'/'`
mydir=`pwd`
mysdb=`echo $mydir | cut -f 7 -d'/'`
if [ $mysrc != $mysdb ]
then
  echo "Your current directory is in $mysdb :\n $mydir"
  echo "Your installation directory is in $mysrc :\n $G4INSTALL"
  echo "On balance, it is more likely you should source setup.csh"
  echo "   in the correct directory tree"
  echo "   than you wish to work across the prod/dev1/dev2 structures"
  echo "\nupdate.sh is stopping\n"
  exit 24
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
${G4STTDIR}/bin/updt.sh
