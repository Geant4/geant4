#! /bin/sh
# $Id: update.sh,v 1.11 2000-06-15 16:52:52 stesting Exp $
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
