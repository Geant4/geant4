#! /bin/sh
# $Id: update.sh,v 1.4 1999-07-28 10:37:21 stesting Exp $
# Edit stt-prod.sdb or stt-dev.sdb and execute.
# Usage: update.sh [-n] [-d]

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
 mkdir $dir
 echo "$dir created."
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
  -d) DIRECTORIES=-d; export DIRECTORIES;;
esac
done

echo

if [ X$NOTHING = X ]
then
  echo ACTUALLY UPDATING
else
  echo "DOING NOTHING (cvs -n)"
fi

if [ X$DIRECTORIES = X ]
then
  echo NOT GETTING NEW DIRECTORIES
else
  echo GETTING NEW DIRECTORIES
fi

cd $G4INSTALL/..
geant4/tests/tools/bin/updt.sh  < geant4/tests/stt-$REF.sdb \
  > $dir/update.log 2>&1
