#----------------------------------------------------------------
# Last update: 30-Jul-2008
#
# CMS hadronic barrel : 3.7 mm Sci, 5 cm Brass.
#                       I am using "Copper" instead of Brass.
#
#  For  Cu :         X0          Lambda     Lambda/X0
#                   1.4353 cm    15.056 cm    10.49 
#  so:
#                    5 cm = 3.48 X0 = 0.3321 lambda
#
# The closest configuration with an overall thickness of about 
# the usual 10 lambda (more precisely, 9.9628 lambda) is the 
# following:  30  active layers (also  30  readout layers).
# The total amount of Scintillator is 11 cm, corresponding 
# (X0 = 42.4 cm, lambda = 79.36 cm) to 0.26 X0 and 0.14 lambda.
#
#----------------------------------------------------------------
#
/random/setSavingFlag 1
#
/run/verbose 1 
/event/verbose 0 
/tracking/verbose 0
#
#/gun/particle e-
/gun/particle pi-
#/gun/particle proton
#
/gun/energy 100 GeV
#
/mydet/absorberMaterial Copper
/mydet/activeMaterial Scintillator
/mydet/isCalHomogeneous 0
#
/mydet/isUnitInLambda 1
/mydet/absorberTotalLength 9.9628
/mydet/calorimeterRadius 4.98140
/mydet/activeLayerNumber 30
/mydet/activeLayerSize 3.7
/mydet/isRadiusUnitInLambda 1
/mydet/radiusBinSize 0.1
/mydet/radiusBinNumber 10
#
###/mydet/setField 4.0 tesla
#
/mydet/update
#
/run/beamOn 5000 
#
