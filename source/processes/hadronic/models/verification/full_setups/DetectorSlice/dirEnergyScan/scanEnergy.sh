#------------------------------------------------------------------------
# Last update: 31-Jul-2008
#
# This shell script run an energy scan, in the beam energy range
# 1 - 30 GeV, for 1 Physics List, 1 beam particle, and 3 types of
# simplified combined calorimeter (cms, atlasBarrel, atlasHec).
#
# Because of the different Birks quenching in scintillator and LAr,
# it is possible to specify a different name of the Physics List 
# for each of the 3 combined calorimeters, for instance:
#   PHYSICS_LIST_CMS=QGSP_BERT-emBirksSci-hadBirksSci
#   PHYSICS_LIST_ATLASbar=QGSP_BERT-emBirksLAr-hadBirksSci
#   PHYSICS_LIST_ATLASend=QGSP_BERT-emBirksLAr-hadBirksLAr
#
# See "***LOOKHERE***" below for the available choices.
#
#------------------------------------------------------------------------
#
#*****************************LOOKHERE***
PHYSICS_LIST_CMS=QGSP_BERT
PHYSICS_LIST_ATLASbar=QGSP_BERT
PHYSICS_LIST_ATLASend=QGSP_BERT
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
mainDetectorSlice-$PHYSICS_LIST_CMS combinedCMS.g4 > output.log-combinedCMS-$PARTICLE-1GeV-$G4_RELEASE-$PHYSICS_LIST_CMS 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASbar combinedATLASbarrel.g4 > output.log-combinedATLASbarrel-$PARTICLE-1GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASbar 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASend combinedATLASendcap.g4 > output.log-combinedATLASendcap-$PARTICLE-1GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASend 2>&1
cd ..
#
date
echo " === 2 GeV === "
cd dir2GeV
mainDetectorSlice-$PHYSICS_LIST_CMS combinedCMS.g4 > output.log-combinedCMS-$PARTICLE-2GeV-$G4_RELEASE-$PHYSICS_LIST_CMS 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASbar combinedATLASbarrel.g4 > output.log-combinedATLASbarrel-$PARTICLE-2GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASbar 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASend combinedATLASendcap.g4 > output.log-combinedATLASendcap-$PARTICLE-2GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASend 2>&1
cd ..
#
date
echo " === 3 GeV === "
cd dir3GeV
mainDetectorSlice-$PHYSICS_LIST_CMS combinedCMS.g4 > output.log-combinedCMS-$PARTICLE-3GeV-$G4_RELEASE-$PHYSICS_LIST_CMS 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASbar combinedATLASbarrel.g4 > output.log-combinedATLASbarrel-$PARTICLE-3GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASbar 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASend combinedATLASendcap.g4 > output.log-combinedATLASendcap-$PARTICLE-3GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASend 2>&1
cd ..
#
date
echo " === 4 GeV === "
cd dir4GeV
mainDetectorSlice-$PHYSICS_LIST_CMS combinedCMS.g4 > output.log-combinedCMS-$PARTICLE-4GeV-$G4_RELEASE-$PHYSICS_LIST_CMS 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASbar combinedATLASbarrel.g4 > output.log-combinedATLASbarrel-$PARTICLE-4GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASbar 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASend combinedATLASendcap.g4 > output.log-combinedATLASendcap-$PARTICLE-4GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASend 2>&1
cd ..
#
date
echo " === 5 GeV === "
cd dir5GeV
mainDetectorSlice-$PHYSICS_LIST_CMS combinedCMS.g4 > output.log-combinedCMS-$PARTICLE-5GeV-$G4_RELEASE-$PHYSICS_LIST_CMS 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASbar combinedATLASbarrel.g4 > output.log-combinedATLASbarrel-$PARTICLE-5GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASbar 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASend combinedATLASendcap.g4 > output.log-combinedATLASendcap-$PARTICLE-5GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASend 2>&1
cd ..
#
date
echo " === 6 GeV === "
cd dir6GeV
mainDetectorSlice-$PHYSICS_LIST_CMS combinedCMS.g4 > output.log-combinedCMS-$PARTICLE-6GeV-$G4_RELEASE-$PHYSICS_LIST_CMS 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASbar combinedATLASbarrel.g4 > output.log-combinedATLASbarrel-$PARTICLE-6GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASbar 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASend combinedATLASendcap.g4 > output.log-combinedATLASendcap-$PARTICLE-6GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASend 2>&1
cd ..
#
date
echo " === 7 GeV === "
cd dir7GeV
mainDetectorSlice-$PHYSICS_LIST_CMS combinedCMS.g4 > output.log-combinedCMS-$PARTICLE-7GeV-$G4_RELEASE-$PHYSICS_LIST_CMS 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASbar combinedATLASbarrel.g4 > output.log-combinedATLASbarrel-$PARTICLE-7GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASbar 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASend combinedATLASendcap.g4 > output.log-combinedATLASendcap-$PARTICLE-7GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASend 2>&1
cd ..
#
date
echo " === 8 GeV === "
cd dir8GeV
mainDetectorSlice-$PHYSICS_LIST_CMS combinedCMS.g4 > output.log-combinedCMS-$PARTICLE-8GeV-$G4_RELEASE-$PHYSICS_LIST_CMS 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASbar combinedATLASbarrel.g4 > output.log-combinedATLASbarrel-$PARTICLE-8GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASbar 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASend combinedATLASendcap.g4 > output.log-combinedATLASendcap-$PARTICLE-8GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASend 2>&1
cd ..
#
date
echo " === 9 GeV === "
cd dir9GeV
mainDetectorSlice-$PHYSICS_LIST_CMS combinedCMS.g4 > output.log-combinedCMS-$PARTICLE-9GeV-$G4_RELEASE-$PHYSICS_LIST_CMS 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASbar combinedATLASbarrel.g4 > output.log-combinedATLASbarrel-$PARTICLE-9GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASbar 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASend combinedATLASendcap.g4 > output.log-combinedATLASendcap-$PARTICLE-9GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASend 2>&1
cd ..
#
date
echo " === 10 GeV === "
cd dir10GeV
mainDetectorSlice-$PHYSICS_LIST_CMS combinedCMS.g4 > output.log-combinedCMS-$PARTICLE-10GeV-$G4_RELEASE-$PHYSICS_LIST_CMS 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASbar combinedATLASbarrel.g4 > output.log-combinedATLASbarrel-$PARTICLE-10GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASbar 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASend combinedATLASendcap.g4 > output.log-combinedATLASendcap-$PARTICLE-10GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASend 2>&1
cd ..
#
date
echo " === 11 GeV === "
cd dir11GeV
mainDetectorSlice-$PHYSICS_LIST_CMS combinedCMS.g4 > output.log-combinedCMS-$PARTICLE-11GeV-$G4_RELEASE-$PHYSICS_LIST_CMS 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASbar combinedATLASbarrel.g4 > output.log-combinedATLASbarrel-$PARTICLE-11GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASbar 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASend combinedATLASendcap.g4 > output.log-combinedATLASendcap-$PARTICLE-11GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASend 2>&1
cd ..
#
date
echo " === 12 GeV === "
cd dir12GeV
mainDetectorSlice-$PHYSICS_LIST_CMS combinedCMS.g4 > output.log-combinedCMS-$PARTICLE-12GeV-$G4_RELEASE-$PHYSICS_LIST_CMS 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASbar combinedATLASbarrel.g4 > output.log-combinedATLASbarrel-$PARTICLE-12GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASbar 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASend combinedATLASendcap.g4 > output.log-combinedATLASendcap-$PARTICLE-12GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASend 2>&1
cd ..
#
date
echo " === 13 GeV === "
cd dir13GeV
mainDetectorSlice-$PHYSICS_LIST_CMS combinedCMS.g4 > output.log-combinedCMS-$PARTICLE-13GeV-$G4_RELEASE-$PHYSICS_LIST_CMS 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASbar combinedATLASbarrel.g4 > output.log-combinedATLASbarrel-$PARTICLE-13GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASbar 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASend combinedATLASendcap.g4 > output.log-combinedATLASendcap-$PARTICLE-13GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASend 2>&1
cd ..
#
date
echo " === 14 GeV === "
cd dir14GeV
mainDetectorSlice-$PHYSICS_LIST_CMS combinedCMS.g4 > output.log-combinedCMS-$PARTICLE-14GeV-$G4_RELEASE-$PHYSICS_LIST_CMS 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASbar combinedATLASbarrel.g4 > output.log-combinedATLASbarrel-$PARTICLE-14GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASbar 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASend combinedATLASendcap.g4 > output.log-combinedATLASendcap-$PARTICLE-14GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASend 2>&1
cd ..
#
date
echo " === 15 GeV === "
cd dir15GeV
mainDetectorSlice-$PHYSICS_LIST_CMS combinedCMS.g4 > output.log-combinedCMS-$PARTICLE-15GeV-$G4_RELEASE-$PHYSICS_LIST_CMS 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASbar combinedATLASbarrel.g4 > output.log-combinedATLASbarrel-$PARTICLE-15GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASbar 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASend combinedATLASendcap.g4 > output.log-combinedATLASendcap-$PARTICLE-15GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASend 2>&1
cd ..
#
date
echo " === 16 GeV === "
cd dir16GeV
mainDetectorSlice-$PHYSICS_LIST_CMS combinedCMS.g4 > output.log-combinedCMS-$PARTICLE-16GeV-$G4_RELEASE-$PHYSICS_LIST_CMS 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASbar combinedATLASbarrel.g4 > output.log-combinedATLASbarrel-$PARTICLE-16GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASbar 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASend combinedATLASendcap.g4 > output.log-combinedATLASendcap-$PARTICLE-16GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASend 2>&1
cd ..
#
date
echo " === 17 GeV === "
cd dir17GeV
mainDetectorSlice-$PHYSICS_LIST_CMS combinedCMS.g4 > output.log-combinedCMS-$PARTICLE-17GeV-$G4_RELEASE-$PHYSICS_LIST_CMS 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASbar combinedATLASbarrel.g4 > output.log-combinedATLASbarrel-$PARTICLE-17GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASbar 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASend combinedATLASendcap.g4 > output.log-combinedATLASendcap-$PARTICLE-17GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASend 2>&1
cd ..
#
date
echo " === 18 GeV === "
cd dir18GeV
mainDetectorSlice-$PHYSICS_LIST_CMS combinedCMS.g4 > output.log-combinedCMS-$PARTICLE-18GeV-$G4_RELEASE-$PHYSICS_LIST_CMS 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASbar combinedATLASbarrel.g4 > output.log-combinedATLASbarrel-$PARTICLE-18GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASbar 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASend combinedATLASendcap.g4 > output.log-combinedATLASendcap-$PARTICLE-18GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASend 2>&1
cd ..
#
date
echo " === 19 GeV === "
cd dir19GeV
mainDetectorSlice-$PHYSICS_LIST_CMS combinedCMS.g4 > output.log-combinedCMS-$PARTICLE-19GeV-$G4_RELEASE-$PHYSICS_LIST_CMS 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASbar combinedATLASbarrel.g4 > output.log-combinedATLASbarrel-$PARTICLE-19GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASbar 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASend combinedATLASendcap.g4 > output.log-combinedATLASendcap-$PARTICLE-19GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASend 2>&1
cd ..
#
date
echo " === 20 GeV === "
cd dir20GeV
mainDetectorSlice-$PHYSICS_LIST_CMS combinedCMS.g4 > output.log-combinedCMS-$PARTICLE-20GeV-$G4_RELEASE-$PHYSICS_LIST_CMS 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASbar combinedATLASbarrel.g4 > output.log-combinedATLASbarrel-$PARTICLE-20GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASbar 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASend combinedATLASendcap.g4 > output.log-combinedATLASendcap-$PARTICLE-20GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASend 2>&1
cd ..
#
date
echo " === 21 GeV === "
cd dir21GeV
mainDetectorSlice-$PHYSICS_LIST_CMS combinedCMS.g4 > output.log-combinedCMS-$PARTICLE-21GeV-$G4_RELEASE-$PHYSICS_LIST_CMS 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASbar combinedATLASbarrel.g4 > output.log-combinedATLASbarrel-$PARTICLE-21GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASbar 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASend combinedATLASendcap.g4 > output.log-combinedATLASendcap-$PARTICLE-21GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASend 2>&1
cd ..
#
date
echo " === 22 GeV === "
cd dir22GeV
mainDetectorSlice-$PHYSICS_LIST_CMS combinedCMS.g4 > output.log-combinedCMS-$PARTICLE-22GeV-$G4_RELEASE-$PHYSICS_LIST_CMS 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASbar combinedATLASbarrel.g4 > output.log-combinedATLASbarrel-$PARTICLE-22GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASbar 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASend combinedATLASendcap.g4 > output.log-combinedATLASendcap-$PARTICLE-22GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASend 2>&1
cd ..
#
date
echo " === 23 GeV === "
cd dir23GeV
mainDetectorSlice-$PHYSICS_LIST_CMS combinedCMS.g4 > output.log-combinedCMS-$PARTICLE-23GeV-$G4_RELEASE-$PHYSICS_LIST_CMS 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASbar combinedATLASbarrel.g4 > output.log-combinedATLASbarrel-$PARTICLE-23GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASbar 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASend combinedATLASendcap.g4 > output.log-combinedATLASendcap-$PARTICLE-23GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASend 2>&1
cd ..
#
date
echo " === 24 GeV === "
cd dir24GeV
mainDetectorSlice-$PHYSICS_LIST_CMS combinedCMS.g4 > output.log-combinedCMS-$PARTICLE-24GeV-$G4_RELEASE-$PHYSICS_LIST_CMS 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASbar combinedATLASbarrel.g4 > output.log-combinedATLASbarrel-$PARTICLE-24GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASbar 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASend combinedATLASendcap.g4 > output.log-combinedATLASendcap-$PARTICLE-24GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASend 2>&1
cd ..
#
date
echo " === 25 GeV === "
cd dir25GeV
mainDetectorSlice-$PHYSICS_LIST_CMS combinedCMS.g4 > output.log-combinedCMS-$PARTICLE-25GeV-$G4_RELEASE-$PHYSICS_LIST_CMS 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASbar combinedATLASbarrel.g4 > output.log-combinedATLASbarrel-$PARTICLE-25GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASbar 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASend combinedATLASendcap.g4 > output.log-combinedATLASendcap-$PARTICLE-25GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASend 2>&1
cd ..
#
date
echo " === 26 GeV === "
cd dir26GeV
mainDetectorSlice-$PHYSICS_LIST_CMS combinedCMS.g4 > output.log-combinedCMS-$PARTICLE-26GeV-$G4_RELEASE-$PHYSICS_LIST_CMS 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASbar combinedATLASbarrel.g4 > output.log-combinedATLASbarrel-$PARTICLE-26GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASbar 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASend combinedATLASendcap.g4 > output.log-combinedATLASendcap-$PARTICLE-26GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASend 2>&1
cd ..
#
date
echo " === 27 GeV === "
cd dir27GeV
mainDetectorSlice-$PHYSICS_LIST_CMS combinedCMS.g4 > output.log-combinedCMS-$PARTICLE-27GeV-$G4_RELEASE-$PHYSICS_LIST_CMS 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASbar combinedATLASbarrel.g4 > output.log-combinedATLASbarrel-$PARTICLE-27GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASbar 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASend combinedATLASendcap.g4 > output.log-combinedATLASendcap-$PARTICLE-27GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASend 2>&1
cd ..
#
date
echo " === 28 GeV === "
cd dir28GeV
mainDetectorSlice-$PHYSICS_LIST_CMS combinedCMS.g4 > output.log-combinedCMS-$PARTICLE-28GeV-$G4_RELEASE-$PHYSICS_LIST_CMS 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASbar combinedATLASbarrel.g4 > output.log-combinedATLASbarrel-$PARTICLE-28GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASbar 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASend combinedATLASendcap.g4 > output.log-combinedATLASendcap-$PARTICLE-28GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASend 2>&1
cd ..
#
date
echo " === 29 GeV === "
cd dir29GeV
mainDetectorSlice-$PHYSICS_LIST_CMS combinedCMS.g4 > output.log-combinedCMS-$PARTICLE-29GeV-$G4_RELEASE-$PHYSICS_LIST_CMS 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASbar combinedATLASbarrel.g4 > output.log-combinedATLASbarrel-$PARTICLE-29GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASbar 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASend combinedATLASendcap.g4 > output.log-combinedATLASendcap-$PARTICLE-29GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASend 2>&1
cd ..
#
date
echo " === 30 GeV === "
cd dir30GeV
mainDetectorSlice-$PHYSICS_LIST_CMS combinedCMS.g4 > output.log-combinedCMS-$PARTICLE-30GeV-$G4_RELEASE-$PHYSICS_LIST_CMS 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASbar combinedATLASbarrel.g4 > output.log-combinedATLASbarrel-$PARTICLE-30GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASbar 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASend combinedATLASendcap.g4 > output.log-combinedATLASendcap-$PARTICLE-30GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASend 2>&1
cd ..
#
echo " "
echo " === END === "
echo " "
