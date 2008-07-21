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
echo " === 1 GeV === "
cd dir1GeV
mainStatAccepTest-$PHYSICS_LIST cmsHCAL.g4 > output.log-cms-$PARTICLE-1GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST atlasHEC.g4 > output.log-atlas-$PARTICLE-1GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST pbwo4.g4 > output.log-pbwo4-$PARTICLE-1GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST tileCal.g4 > output.log-tile-$PARTICLE-1GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 2 GeV === "
cd dir2GeV
mainStatAccepTest-$PHYSICS_LIST cmsHCAL.g4 > output.log-cms-$PARTICLE-2GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST atlasHEC.g4 > output.log-atlas-$PARTICLE-2GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST pbwo4.g4 > output.log-pbwo4-$PARTICLE-2GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST tileCal.g4 > output.log-tile-$PARTICLE-2GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 3 GeV === "
cd dir3GeV
mainStatAccepTest-$PHYSICS_LIST cmsHCAL.g4 > output.log-cms-$PARTICLE-3GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST atlasHEC.g4 > output.log-atlas-$PARTICLE-3GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST pbwo4.g4 > output.log-pbwo4-$PARTICLE-3GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST tileCal.g4 > output.log-tile-$PARTICLE-3GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 4 GeV === "
cd dir4GeV
mainStatAccepTest-$PHYSICS_LIST cmsHCAL.g4 > output.log-cms-$PARTICLE-4GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST atlasHEC.g4 > output.log-atlas-$PARTICLE-4GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST pbwo4.g4 > output.log-pbwo4-$PARTICLE-4GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST tileCal.g4 > output.log-tile-$PARTICLE-4GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 5 GeV === "
cd dir5GeV
mainStatAccepTest-$PHYSICS_LIST cmsHCAL.g4 > output.log-cms-$PARTICLE-5GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST atlasHEC.g4 > output.log-atlas-$PARTICLE-5GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST pbwo4.g4 > output.log-pbwo4-$PARTICLE-5GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST tileCal.g4 > output.log-tile-$PARTICLE-5GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 6 GeV === "
cd dir6GeV
mainStatAccepTest-$PHYSICS_LIST cmsHCAL.g4 > output.log-cms-$PARTICLE-6GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST atlasHEC.g4 > output.log-atlas-$PARTICLE-6GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST pbwo4.g4 > output.log-pbwo4-$PARTICLE-6GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST tileCal.g4 > output.log-tile-$PARTICLE-6GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 7 GeV === "
cd dir7GeV
mainStatAccepTest-$PHYSICS_LIST cmsHCAL.g4 > output.log-cms-$PARTICLE-7GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST atlasHEC.g4 > output.log-atlas-$PARTICLE-7GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST pbwo4.g4 > output.log-pbwo4-$PARTICLE-7GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST tileCal.g4 > output.log-tile-$PARTICLE-7GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 8 GeV === "
cd dir8GeV
mainStatAccepTest-$PHYSICS_LIST cmsHCAL.g4 > output.log-cms-$PARTICLE-8GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST atlasHEC.g4 > output.log-atlas-$PARTICLE-8GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST pbwo4.g4 > output.log-pbwo4-$PARTICLE-8GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST tileCal.g4 > output.log-tile-$PARTICLE-8GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 9 GeV === "
cd dir9GeV
mainStatAccepTest-$PHYSICS_LIST cmsHCAL.g4 > output.log-cms-$PARTICLE-9GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST atlasHEC.g4 > output.log-atlas-$PARTICLE-9GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST pbwo4.g4 > output.log-pbwo4-$PARTICLE-9GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST tileCal.g4 > output.log-tile-$PARTICLE-9GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 10 GeV === "
cd dir10GeV
mainStatAccepTest-$PHYSICS_LIST cmsHCAL.g4 > output.log-cms-$PARTICLE-10GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST atlasHEC.g4 > output.log-atlas-$PARTICLE-10GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST pbwo4.g4 > output.log-pbwo4-$PARTICLE-10GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST tileCal.g4 > output.log-tile-$PARTICLE-10GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 11 GeV === "
cd dir11GeV
mainStatAccepTest-$PHYSICS_LIST cmsHCAL.g4 > output.log-cms-$PARTICLE-11GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST atlasHEC.g4 > output.log-atlas-$PARTICLE-11GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST pbwo4.g4 > output.log-pbwo4-$PARTICLE-11GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST tileCal.g4 > output.log-tile-$PARTICLE-11GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 12 GeV === "
cd dir12GeV
mainStatAccepTest-$PHYSICS_LIST cmsHCAL.g4 > output.log-cms-$PARTICLE-12GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST atlasHEC.g4 > output.log-atlas-$PARTICLE-12GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST pbwo4.g4 > output.log-pbwo4-$PARTICLE-12GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST tileCal.g4 > output.log-tile-$PARTICLE-12GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 13 GeV === "
cd dir13GeV
mainStatAccepTest-$PHYSICS_LIST cmsHCAL.g4 > output.log-cms-$PARTICLE-13GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST atlasHEC.g4 > output.log-atlas-$PARTICLE-13GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST pbwo4.g4 > output.log-pbwo4-$PARTICLE-13GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST tileCal.g4 > output.log-tile-$PARTICLE-13GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 14 GeV === "
cd dir14GeV
mainStatAccepTest-$PHYSICS_LIST cmsHCAL.g4 > output.log-cms-$PARTICLE-14GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST atlasHEC.g4 > output.log-atlas-$PARTICLE-14GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST pbwo4.g4 > output.log-pbwo4-$PARTICLE-14GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST tileCal.g4 > output.log-tile-$PARTICLE-14GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 15 GeV === "
cd dir15GeV
mainStatAccepTest-$PHYSICS_LIST cmsHCAL.g4 > output.log-cms-$PARTICLE-15GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST atlasHEC.g4 > output.log-atlas-$PARTICLE-15GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST pbwo4.g4 > output.log-pbwo4-$PARTICLE-15GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST tileCal.g4 > output.log-tile-$PARTICLE-15GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 16 GeV === "
cd dir16GeV
mainStatAccepTest-$PHYSICS_LIST cmsHCAL.g4 > output.log-cms-$PARTICLE-16GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST atlasHEC.g4 > output.log-atlas-$PARTICLE-16GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST pbwo4.g4 > output.log-pbwo4-$PARTICLE-16GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST tileCal.g4 > output.log-tile-$PARTICLE-16GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 17 GeV === "
cd dir17GeV
mainStatAccepTest-$PHYSICS_LIST cmsHCAL.g4 > output.log-cms-$PARTICLE-17GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST atlasHEC.g4 > output.log-atlas-$PARTICLE-17GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST pbwo4.g4 > output.log-pbwo4-$PARTICLE-17GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST tileCal.g4 > output.log-tile-$PARTICLE-17GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 18 GeV === "
cd dir18GeV
mainStatAccepTest-$PHYSICS_LIST cmsHCAL.g4 > output.log-cms-$PARTICLE-18GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST atlasHEC.g4 > output.log-atlas-$PARTICLE-18GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST pbwo4.g4 > output.log-pbwo4-$PARTICLE-18GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST tileCal.g4 > output.log-tile-$PARTICLE-18GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 19 GeV === "
cd dir19GeV
mainStatAccepTest-$PHYSICS_LIST cmsHCAL.g4 > output.log-cms-$PARTICLE-19GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST atlasHEC.g4 > output.log-atlas-$PARTICLE-19GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST pbwo4.g4 > output.log-pbwo4-$PARTICLE-19GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST tileCal.g4 > output.log-tile-$PARTICLE-19GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 20 GeV === "
cd dir20GeV
mainStatAccepTest-$PHYSICS_LIST cmsHCAL.g4 > output.log-cms-$PARTICLE-20GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST atlasHEC.g4 > output.log-atlas-$PARTICLE-20GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST pbwo4.g4 > output.log-pbwo4-$PARTICLE-20GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST tileCal.g4 > output.log-tile-$PARTICLE-20GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 21 GeV === "
cd dir21GeV
mainStatAccepTest-$PHYSICS_LIST cmsHCAL.g4 > output.log-cms-$PARTICLE-21GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST atlasHEC.g4 > output.log-atlas-$PARTICLE-21GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST pbwo4.g4 > output.log-pbwo4-$PARTICLE-21GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST tileCal.g4 > output.log-tile-$PARTICLE-21GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 22 GeV === "
cd dir22GeV
mainStatAccepTest-$PHYSICS_LIST cmsHCAL.g4 > output.log-cms-$PARTICLE-22GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST atlasHEC.g4 > output.log-atlas-$PARTICLE-22GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST pbwo4.g4 > output.log-pbwo4-$PARTICLE-22GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST tileCal.g4 > output.log-tile-$PARTICLE-22GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 23 GeV === "
cd dir23GeV
mainStatAccepTest-$PHYSICS_LIST cmsHCAL.g4 > output.log-cms-$PARTICLE-23GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST atlasHEC.g4 > output.log-atlas-$PARTICLE-23GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST pbwo4.g4 > output.log-pbwo4-$PARTICLE-23GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST tileCal.g4 > output.log-tile-$PARTICLE-23GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 24 GeV === "
cd dir24GeV
mainStatAccepTest-$PHYSICS_LIST cmsHCAL.g4 > output.log-cms-$PARTICLE-24GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST atlasHEC.g4 > output.log-atlas-$PARTICLE-24GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST pbwo4.g4 > output.log-pbwo4-$PARTICLE-24GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST tileCal.g4 > output.log-tile-$PARTICLE-24GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 25 GeV === "
cd dir25GeV
mainStatAccepTest-$PHYSICS_LIST cmsHCAL.g4 > output.log-cms-$PARTICLE-25GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST atlasHEC.g4 > output.log-atlas-$PARTICLE-25GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST pbwo4.g4 > output.log-pbwo4-$PARTICLE-25GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST tileCal.g4 > output.log-tile-$PARTICLE-25GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 26 GeV === "
cd dir26GeV
mainStatAccepTest-$PHYSICS_LIST cmsHCAL.g4 > output.log-cms-$PARTICLE-26GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST atlasHEC.g4 > output.log-atlas-$PARTICLE-26GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST pbwo4.g4 > output.log-pbwo4-$PARTICLE-26GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST tileCal.g4 > output.log-tile-$PARTICLE-26GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 27 GeV === "
cd dir27GeV
mainStatAccepTest-$PHYSICS_LIST cmsHCAL.g4 > output.log-cms-$PARTICLE-27GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST atlasHEC.g4 > output.log-atlas-$PARTICLE-27GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST pbwo4.g4 > output.log-pbwo4-$PARTICLE-27GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST tileCal.g4 > output.log-tile-$PARTICLE-27GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 28 GeV === "
cd dir28GeV
mainStatAccepTest-$PHYSICS_LIST cmsHCAL.g4 > output.log-cms-$PARTICLE-28GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST atlasHEC.g4 > output.log-atlas-$PARTICLE-28GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST pbwo4.g4 > output.log-pbwo4-$PARTICLE-28GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST tileCal.g4 > output.log-tile-$PARTICLE-28GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 29 GeV === "
cd dir29GeV
mainStatAccepTest-$PHYSICS_LIST cmsHCAL.g4 > output.log-cms-$PARTICLE-29GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST atlasHEC.g4 > output.log-atlas-$PARTICLE-29GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST pbwo4.g4 > output.log-pbwo4-$PARTICLE-29GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST tileCal.g4 > output.log-tile-$PARTICLE-29GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
date
echo " === 30 GeV === "
cd dir30GeV
mainStatAccepTest-$PHYSICS_LIST cmsHCAL.g4 > output.log-cms-$PARTICLE-30GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST atlasHEC.g4 > output.log-atlas-$PARTICLE-30GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST pbwo4.g4 > output.log-pbwo4-$PARTICLE-30GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
mainStatAccepTest-$PHYSICS_LIST tileCal.g4 > output.log-tile-$PARTICLE-30GeV-$G4_RELEASE-$PHYSICS_LIST 2>&1
cd ..
#
echo " "
echo " === END === "
echo " "
