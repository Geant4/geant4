#------------------------------------------------------------------------
# Last update: 21-Aug-2008
#
# This simple shell script run a energy scan, for few beam energy
# points above 30 GeV, for 1 Physics List, 1 beam particle, and 
# 1 type of simplified calorimeter.
#
# See "***LOOKHERE***" below for the available choices.
#
# This script is similar to  extraScanEnergy.sh : the difference is 
# that the latter runs 4 simplified calorimeters (cms, atlas, pbwo4, tile),
# whereas  extraScanEnergy1.sh  runs only 1 simplified calorimeter.
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
#
###CALO=cmsHCAL   ; sCALO=cms
###CALO=atlasHEC  ; sCALO=atlas
###CALO=pbwo4     ; sCALO=pbwo4
###CALO=tileCal   ; sCALO=tile
CALO=atlasTILE ; sCALO=atile
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
mainStatAccepTest-$PHYSICS_LIST $CALO.g4 > output.log-$sCALO-$PARTICLE-50GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 100 GeV === "
cd dir100GeV
mainStatAccepTest-$PHYSICS_LIST $CALO.g4 > output.log-$sCALO-$PARTICLE-100GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
###date
###echo " === 180 GeV === "
###cd dir180GeV
###mainStatAccepTest-$PHYSICS_LIST $CALO.g4 > output.log-$sCALO-$PARTICLE-180GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
###cd ..
#
date
echo " === 200 GeV === "
cd dir200GeV
mainStatAccepTest-$PHYSICS_LIST $CALO.g4 > output.log-$sCALO-$PARTICLE-200GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
###date
###echo " === 300 GeV === "
###cd dir300GeV
###mainStatAccepTest-$PHYSICS_LIST $CALO.g4 > output.log-$sCALO-$PARTICLE-300GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
###cd ..
#
###date
###echo " === 350 GeV === "
###cd dir350GeV
###mainStatAccepTest-$PHYSICS_LIST $CALO.g4 > output.log-$sCALO-$PARTICLE-350GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
###cd ..
#
date
#
echo " "
echo " === END === "
echo " "
