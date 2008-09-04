#------------------------------------------------------------------------
# Last update: 04-Sep-2008
#
# This simple shell script run a energy scan, in the beam energy range
# 1 - 30 GeV, and beyond, for 1 Physics List, 1 beam particle, 
# and 1 type of simplified calorimeter.
#
# See "***LOOKHERE***" below for the available choices.
#
# This script is similar to  scanEnergy.sh : the difference is that
# the latter runs 4 simplified calorimeters (cms, atlas, pbwo4, tile),
# whereas  scanEnergy1.sh  runs only 1 simplified calorimeter.
# Notice that  scanEnergy1.sh  can run  atlasTILE.g4 , whereas
# scanEnergy.sh  cannot. 
# Furthermore, it can run also beam energies above 30 GeV.
#
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
echo " === 1 GeV === "
cd dir1GeV
mainStatAccepTest-$PHYSICS_LIST $CALO.g4 > output.log-$sCALO-$PARTICLE-1GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 2 GeV === "
cd dir2GeV
mainStatAccepTest-$PHYSICS_LIST $CALO.g4 > output.log-$sCALO-$PARTICLE-2GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 3 GeV === "
cd dir3GeV
mainStatAccepTest-$PHYSICS_LIST $CALO.g4 > output.log-$sCALO-$PARTICLE-3GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 4 GeV === "
cd dir4GeV
mainStatAccepTest-$PHYSICS_LIST $CALO.g4 > output.log-$sCALO-$PARTICLE-4GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 5 GeV === "
cd dir5GeV
mainStatAccepTest-$PHYSICS_LIST $CALO.g4 > output.log-$sCALO-$PARTICLE-5GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 6 GeV === "
cd dir6GeV
mainStatAccepTest-$PHYSICS_LIST $CALO.g4 > output.log-$sCALO-$PARTICLE-6GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 7 GeV === "
cd dir7GeV
mainStatAccepTest-$PHYSICS_LIST $CALO.g4 > output.log-$sCALO-$PARTICLE-7GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 8 GeV === "
cd dir8GeV
mainStatAccepTest-$PHYSICS_LIST $CALO.g4 > output.log-$sCALO-$PARTICLE-8GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 9 GeV === "
cd dir9GeV
mainStatAccepTest-$PHYSICS_LIST $CALO.g4 > output.log-$sCALO-$PARTICLE-9GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 10 GeV === "
cd dir10GeV
mainStatAccepTest-$PHYSICS_LIST $CALO.g4 > output.log-$sCALO-$PARTICLE-10GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 11 GeV === "
cd dir11GeV
mainStatAccepTest-$PHYSICS_LIST $CALO.g4 > output.log-$sCALO-$PARTICLE-11GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 12 GeV === "
cd dir12GeV
mainStatAccepTest-$PHYSICS_LIST $CALO.g4 > output.log-$sCALO-$PARTICLE-12GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 13 GeV === "
cd dir13GeV
mainStatAccepTest-$PHYSICS_LIST $CALO.g4 > output.log-$sCALO-$PARTICLE-13GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 14 GeV === "
cd dir14GeV
mainStatAccepTest-$PHYSICS_LIST $CALO.g4 > output.log-$sCALO-$PARTICLE-14GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 15 GeV === "
cd dir15GeV
mainStatAccepTest-$PHYSICS_LIST $CALO.g4 > output.log-$sCALO-$PARTICLE-15GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 16 GeV === "
cd dir16GeV
mainStatAccepTest-$PHYSICS_LIST $CALO.g4 > output.log-$sCALO-$PARTICLE-16GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 17 GeV === "
cd dir17GeV
mainStatAccepTest-$PHYSICS_LIST $CALO.g4 > output.log-$sCALO-$PARTICLE-17GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 18 GeV === "
cd dir18GeV
mainStatAccepTest-$PHYSICS_LIST $CALO.g4 > output.log-$sCALO-$PARTICLE-18GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 19 GeV === "
cd dir19GeV
mainStatAccepTest-$PHYSICS_LIST $CALO.g4 > output.log-$sCALO-$PARTICLE-19GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 20 GeV === "
cd dir20GeV
mainStatAccepTest-$PHYSICS_LIST $CALO.g4 > output.log-$sCALO-$PARTICLE-20GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 21 GeV === "
cd dir21GeV
mainStatAccepTest-$PHYSICS_LIST $CALO.g4 > output.log-$sCALO-$PARTICLE-21GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 22 GeV === "
cd dir22GeV
mainStatAccepTest-$PHYSICS_LIST $CALO.g4 > output.log-$sCALO-$PARTICLE-22GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 23 GeV === "
cd dir23GeV
mainStatAccepTest-$PHYSICS_LIST $CALO.g4 > output.log-$sCALO-$PARTICLE-23GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 24 GeV === "
cd dir24GeV
mainStatAccepTest-$PHYSICS_LIST $CALO.g4 > output.log-$sCALO-$PARTICLE-24GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 25 GeV === "
cd dir25GeV
mainStatAccepTest-$PHYSICS_LIST $CALO.g4 > output.log-$sCALO-$PARTICLE-25GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 26 GeV === "
cd dir26GeV
mainStatAccepTest-$PHYSICS_LIST $CALO.g4 > output.log-$sCALO-$PARTICLE-26GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 27 GeV === "
cd dir27GeV
mainStatAccepTest-$PHYSICS_LIST $CALO.g4 > output.log-$sCALO-$PARTICLE-27GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 28 GeV === "
cd dir28GeV
mainStatAccepTest-$PHYSICS_LIST $CALO.g4 > output.log-$sCALO-$PARTICLE-28GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 29 GeV === "
cd dir29GeV
mainStatAccepTest-$PHYSICS_LIST $CALO.g4 > output.log-$sCALO-$PARTICLE-29GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 30 GeV === "
cd dir30GeV
mainStatAccepTest-$PHYSICS_LIST $CALO.g4 > output.log-$sCALO-$PARTICLE-30GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
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
