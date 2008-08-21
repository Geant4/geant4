#------------------------------------------------------------------------
# Last update: 21-Aug-2008
#
# This simple shell script run a energy scan, for few beam energy
# points above 30 GeV, for 1 Physics List, 1 beam particle, and 
# 4 types of simplified calorimeter (cms, atlas, pbwo4, tile).
#
# See "***LOOKHERE***" below for the available choices.
#
# This script is similar to  extraScanEnergy1.sh : the difference is 
# that the latter runs only 1 simplified calorimeter, whereas  scanEnergy.sh
# runs 4 simplified calorimeters (cms, atlas, pbwo4, tile),
# Notice that  extraScanEnergy1.sh  can run  atlasTILE.g4 , whereas
# extraScanEnergy.sh  cannot.
#------------------------------------------------------------------------
#
#*****************************LOOKHERE***
PHYSICS_LIST=QGSP_BERT
#
G4_RELEASE=9.1.p02
#
###PARTICLE=e
PARTICLE=pi
###PARTICLE=pi+
###PARTICLE=p
#****************************************
#
cd dirEnergyScan_$PARTICLE
#
echo " "
echo " === BEGIN === " $PHYSICS_LIST
echo " "
#
date
echo " === 50 GeV === "
cd dir50GeV
mainStatAccepTest-$PHYSICS_LIST cmsHCAL.g4 > output.log-cms-$PARTICLE-50GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST atlasHEC.g4 > output.log-atlas-$PARTICLE-50GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST pbwo4.g4 > output.log-pbwo4-$PARTICLE-50GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST tileCal.g4 > output.log-tile-$PARTICLE-50GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 100 GeV === "
cd dir100GeV
mainStatAccepTest-$PHYSICS_LIST cmsHCAL.g4 > output.log-cms-$PARTICLE-100GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST atlasHEC.g4 > output.log-atlas-$PARTICLE-100GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST pbwo4.g4 > output.log-pbwo4-$PARTICLE-100GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST tileCal.g4 > output.log-tile-$PARTICLE-100GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
###date
###echo " === 180 GeV === "
###cd dir180GeV
###mainStatAccepTest-$PHYSICS_LIST cmsHCAL.g4 > output.log-cms-$PARTICLE-180GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
###mainStatAccepTest-$PHYSICS_LIST atlasHEC.g4 > output.log-atlas-$PARTICLE-180GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
###mainStatAccepTest-$PHYSICS_LIST pbwo4.g4 > output.log-pbwo4-$PARTICLE-180GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
###mainStatAccepTest-$PHYSICS_LIST tileCal.g4 > output.log-tile-$PARTICLE-180GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
###cd ..
#
date
echo " === 200 GeV === "
cd dir200GeV
mainStatAccepTest-$PHYSICS_LIST cmsHCAL.g4 > output.log-cms-$PARTICLE-200GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST atlasHEC.g4 > output.log-atlas-$PARTICLE-200GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST pbwo4.g4 > output.log-pbwo4-$PARTICLE-200GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST tileCal.g4 > output.log-tile-$PARTICLE-200GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
###date
###echo " === 300 GeV === "
###cd dir300GeV
###mainStatAccepTest-$PHYSICS_LIST cmsHCAL.g4 > output.log-cms-$PARTICLE-300GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
###mainStatAccepTest-$PHYSICS_LIST atlasHEC.g4 > output.log-atlas-$PARTICLE-300GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
###mainStatAccepTest-$PHYSICS_LIST pbwo4.g4 > output.log-pbwo4-$PARTICLE-300GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
###mainStatAccepTest-$PHYSICS_LIST tileCal.g4 > output.log-tile-$PARTICLE-300GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
###cd ..
#
###date
###echo " === 350 GeV === "
###cd dir350GeV
###mainStatAccepTest-$PHYSICS_LIST cmsHCAL.g4 > output.log-cms-$PARTICLE-350GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
###mainStatAccepTest-$PHYSICS_LIST atlasHEC.g4 > output.log-atlas-$PARTICLE-350GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
###mainStatAccepTest-$PHYSICS_LIST pbwo4.g4 > output.log-pbwo4-$PARTICLE-350GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
###mainStatAccepTest-$PHYSICS_LIST tileCal.g4 > output.log-tile-$PARTICLE-350GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
###cd ..
#
date
#
echo " "
echo " === END === "
echo " "
