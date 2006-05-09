#!/bin/csh -f
#trace on
#---------------------------------------------------------------
#
#  Author V.Ivanchenko 26 May 2005
#  Standard EM test for hadronic generator of HARP
#----------------------------------------------------------------

#cd $G4INSTALL/source
#gmake
source ~/bin/pi_setup.csh
cd $G4INSTALL/source/processes/hadronic/models/verification/thin_target/harp
setenv TARGET p_al_13gev
gmake
echo "Start of run for " $TARGET
mkdir $1
cd $1

#$G4MY/harp ../$TARGET/run.mac  >&  res.out
$G4MY/harp ../$TARGET/runMay05.mac >& res.out

cd ../

echo $TARGET " is done!"
