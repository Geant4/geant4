#------------------------------------------------------------------------
# Last update: 21-Aug-2008
#
# This shell script run an energy scan, for some beam energies
# above 30 GeV, for 1 Physics List, 1 beam particle, and 3 types of
# simplified combined calorimeter (cms, atlasBarrel, atlasHec),
# with only the HAD calorimeters (i.e. without their EM ones).
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
PHYSICS_LIST_CMS=QGSP_BERT-hadBirks
PHYSICS_LIST_ATLASbar=QGSP_BERT-hadBirks
PHYSICS_LIST_ATLASend=QGSP_BERT
#
G4_RELEASE=9.2.b01
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
echo " === 50 GeV === "
cd dir50GeV
mainDetectorSlice-$PHYSICS_LIST_CMS combinedCMS.g4 > output.log-combinedCMS_onlyHad-$PARTICLE-50GeV-$G4_RELEASE-$PHYSICS_LIST_CMS 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASbar combinedATLASbarrel.g4 > output.log-combinedATLASbarrel_onlyHad-$PARTICLE-50GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASbar 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASend combinedATLASendcap.g4 > output.log-combinedATLASendcap_onlyHad-$PARTICLE-50GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASend 2>&1
cd ..
#
echo " === 100 GeV === "
cd dir100GeV
mainDetectorSlice-$PHYSICS_LIST_CMS combinedCMS.g4 > output.log-combinedCMS_onlyHad-$PARTICLE-100GeV-$G4_RELEASE-$PHYSICS_LIST_CMS 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASbar combinedATLASbarrel.g4 > output.log-combinedATLASbarrel_onlyHad-$PARTICLE-100GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASbar 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASend combinedATLASendcap.g4 > output.log-combinedATLASendcap_onlyHad-$PARTICLE-100GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASend 2>&1
cd ..
#
###echo " === 180 GeV === "
###cd dir180GeV
###mainDetectorSlice-$PHYSICS_LIST_CMS combinedCMS.g4 > output.log-combinedCMS_onlyHad-$PARTICLE-180GeV-$G4_RELEASE-$PHYSICS_LIST_CMS 2>&1
###mainDetectorSlice-$PHYSICS_LIST_ATLASbar combinedATLASbarrel.g4 > output.log-combinedATLASbarrel_onlyHad-$PARTICLE-180GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASbar 2>&1
###mainDetectorSlice-$PHYSICS_LIST_ATLASend combinedATLASendcap.g4 > output.log-combinedATLASendcap_onlyHad-$PARTICLE-180GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASend 2>&1
###cd ..
#
echo " === 200 GeV === "
cd dir200GeV
mainDetectorSlice-$PHYSICS_LIST_CMS combinedCMS.g4 > output.log-combinedCMS_onlyHad-$PARTICLE-200GeV-$G4_RELEASE-$PHYSICS_LIST_CMS 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASbar combinedATLASbarrel.g4 > output.log-combinedATLASbarrel_onlyHad-$PARTICLE-200GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASbar 2>&1
mainDetectorSlice-$PHYSICS_LIST_ATLASend combinedATLASendcap.g4 > output.log-combinedATLASendcap_onlyHad-$PARTICLE-200GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASend 2>&1
cd ..
#
###echo " === 300 GeV === "
###cd dir300GeV
###mainDetectorSlice-$PHYSICS_LIST_CMS combinedCMS.g4 > output.log-combinedCMS_onlyHad-$PARTICLE-300GeV-$G4_RELEASE-$PHYSICS_LIST_CMS 2>&1
###mainDetectorSlice-$PHYSICS_LIST_ATLASbar combinedATLASbarrel.g4 > output.log-combinedATLASbarrel_onlyHad-$PARTICLE-300GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASbar 2>&1
###mainDetectorSlice-$PHYSICS_LIST_ATLASend combinedATLASendcap.g4 > output.log-combinedATLASendcap_onlyHad-$PARTICLE-300GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASend 2>&1
###cd ..
#
###echo " === 350 GeV === "
###cd dir350GeV
###mainDetectorSlice-$PHYSICS_LIST_CMS combinedCMS.g4 > output.log-combinedCMS_onlyHad-$PARTICLE-350GeV-$G4_RELEASE-$PHYSICS_LIST_CMS 2>&1
###mainDetectorSlice-$PHYSICS_LIST_ATLASbar combinedATLASbarrel.g4 > output.log-combinedATLASbarrel_onlyHad-$PARTICLE-350GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASbar 2>&1
###mainDetectorSlice-$PHYSICS_LIST_ATLASend combinedATLASendcap.g4 > output.log-combinedATLASendcap_onlyHad-$PARTICLE-350GeV-$G4_RELEASE-$PHYSICS_LIST_ATLASend 2>&1
###cd ..
#
echo " "
echo " === END === "
echo " "
