#----------------------------------------------------------------
# Last update: 14-Jun-2006
#
# This configuration, 10 lambda homogeneous PbWO4 calorimeter, 
# with 100 fake layers and 20 readout layers, is meant only 
# for hadronic physics testing purposes (with the same material, 
# but not the geometry, of CMS crystal ECAL).
#
#----------------------------------------------------------------
#
/random/setSavingFlag 1
#
/run/verbose 1 
/event/verbose 0 
/tracking/verbose 0 
#
/gun/particle pi-
#
/gun/energy 100 GeV
#
/mydet/absorberMaterial PbWO4
/mydet/activeMaterial PbWO4
/mydet/isCalHomogeneous 1
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
