#----------------------------------------------------------------
# Last update: 28-Apr-2006
#
# CMS hadronic barrel : 4 mm Sci.
#                       For the absorber thickness, I have found
#			different numbers:
#                         -  5 cm Cu (barrel) and 
#                            8 cm (endcap) 
#                            from an old CMS TRD.
#                         -  6 cm   Cu layers 2-10 and 
#                            6.6 cm Cu layers 11-16 
#                            (layer 1 : steel 7.45 cm thick;
#                             layer 17: steel 8.9  cm thick) 
#                            in Elvira's CMS note.
#                       I assume:  6 cm .
#
#  For  Cu :         X0          Lambda     Lambda/X0
#                   1.4353 cm    15.056 cm    10.49 
#  so:
#                    6 cm = 4.18 X0 = 0.3985 lambda
#
# The closest configuration with an overall thickness of about 
# the usual 10 lambda (more precisely, 9.9625 lambda) is the 
# following:  25  active layers
# which corresponds to  0.3985 lambda  per passive layer.           
# We use  25 readout layers  (which is close to the usual 20). 
# The total amount of Scintillator is 10 cm, corresponding 
# (X0 = 42.4 cm, lambda = 79.36 cm) to 0.24 X0 and 0.13 lambda.
#
# The real setup, as in Elvira's note, is made of layers with
# different thicknesses and even materials, and are in total 17.
# The total thickness corresponds to about  7.2  lambda.
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
#
/gun/energy 100 GeV
#
/mydet/absorberMaterial Copper
/mydet/activeMaterial Scintillator
/mydet/isCalHomogeneous 0
#
/mydet/isUnitInLambda 1
/mydet/absorberTotalLength 9.9625
/mydet/calorimeterRadius 4.98125
/mydet/activeLayerNumber 25
/mydet/activeLayerSize 4.0
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
