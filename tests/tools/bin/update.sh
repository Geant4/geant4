#! /bin/sh
# $Id: update.sh,v 1.1 1999-01-08 16:36:08 gunter Exp $
# Edit stt-ref.sdb and execute.
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

cd $G4INSTALL/..
if [ `pwd | grep ref+` ]
then
  REF=ref+
else
  REF=ref
fi
geant4beta/tests/tools/bin/updt.sh $1 < geant4beta/tests/stt-${REF}.sdb \
  > $dir/update.log 2>&1
