#----------------------------------------------------------------
# Last update: 01-December-2006
#
# This macro file can be used to test Geant4 High Precision (HP)
# treatment of low-energy neutrons (Ekin < 20 MeV), for instance
# with the Physics List QGSP_BERT_HP.
# The absorber is Lead, and you can choose between Liquid Argon
# and Scintillator as active materials.
# The beam particle energy is 3 GeV (similar to TARC), and you
# can choose between proton and pion- as beam particle types.
#
#----------------------------------------------------------------
#
#
/random/setSavingFlag 1
#
/run/verbose 1 
/event/verbose 0 
/tracking/verbose 0 
#
#--- Particle ---
/gun/particle pi-
###/gun/particle proton
#
#
/gun/energy 3 GeV
#
#
/mydet/absorberMaterial Lead
#
#--- Active layer ---
/mydet/activeMaterial LiquidArgon
###/mydet/activeMaterial Scintillator
#
#
/mydet/isCalHomogeneous 0
#
/mydet/isUnitInLambda 1
/mydet/absorberTotalLength 10.0
/mydet/calorimeterRadius 5.0
/mydet/activeLayerNumber 100
/mydet/readoutLayerNumber 20
/mydet/activeLayerSize 4.0
/mydet/isRadiusUnitInLambda 1
/mydet/radiusBinSize 0.1
/mydet/radiusBinNumber 10
#
/mydet/update
#
/run/beamOn 5000
#
