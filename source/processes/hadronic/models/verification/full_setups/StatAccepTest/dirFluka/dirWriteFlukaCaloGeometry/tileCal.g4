#--------------------------------------------------------------------
# Last update: 09-Apr-2006
#
# Both ATLAS TileCal and LHCb hadronic calorimeters are sampling 
# calorimeter  Fe-Sci , but with the scintillator tiles running 
# parallel to the beam, so it is not possible reproduce this 
# configuration with the simplified calorimeters.
# 
# Just to have a reference case for  Fe-Sci , but with the
# scintillator orthogonal to the beam direction, we then use the
# usual setup: 10 lambda calorimeter, with 100 active layers, and 
# 20 readout layers.
#
# For  Fe :         X0          Lambda     Lambda/X0
#                  1.7585 cm    16.760 cm     9.53    
#  so:
#       1 iron layer = 0.1 lambda = 1.676 cm = 0.953 X0
#
#--------------------------------------------------------------------
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
/mydet/absorberMaterial Iron
/mydet/activeMaterial Scintillator
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
