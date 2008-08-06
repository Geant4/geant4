#--------------------------------------------------------------------
# Last update: 28-Apr-2006
#
# ATLAS HEC : HEC1 : Cu (25 mm) LAr (8.5 mm)  0.82 m length
#             HEC2 : Cu (50 mm) LAr (8.5 mm)  0.96 m length
#             4 longitudinal compartments, with a total detector
#             thickness of about  103 X0  and  > 10 lambda.
#
# For  Cu :         X0          Lambda     Lambda/X0
#                  1.4353 cm    15.056 cm    10.49 
# so:
#                  25 mm = 1.74 X0 = 0.1660 lambda
#
# The closest configuration with an overall thickness of about 
# the usual 10 lambda (more precisely, 9.96 lambda) is the 
# following:  60  active layers
# which corresponds to  0.1660 lambda  per passive layer.           
# We use as usual  20 readout layers .
# The total amount of LAr is 51 cm, corresponding 
# (X0 = 10.971 cm, lambda = 65.769 cm) to 4.64 X0 and 0.78 lambda.
#
# The real setup is made of two sets of layers with different 
# thicknesses, and the total thickness seems to me less than
# the claimed > 10 lambda.
# However, a reasonable way to emulate the 4 segmentations
# would be the following: 
#   o  HEC1  I :  1-12 active layers, so  1-4 readout ones;
#   o  HEC1 II : 13-24 active layers, so  5-8 readout ones;
#   o  HEC2  I : 25-42 active layers, so half of  9-15 readout ones; 
#   o  HEC2 II : 43-60 active layers, so half of 16-20 readout ones.
#
# NB) If one wants to use this configuration to get also some 
#     rough information regarding HEC when there is the EMEC 
#     in front of it, the easiest way is the following.
#     The EMEC has about 26 X0, which corresponds roughly to 
#     1 lambda. This corresponds to the first 6 active layers
#     of the HEC, i.e. the first 2 readout layers.
#     Therefore, by looking at the totat energy deposit excluding
#     the first 2 readout layers, one can get a rough estimate of
#     the signal in the HEC in the EMEC+HEC configuration.
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
/mydet/absorberMaterial Copper
/mydet/activeMaterial LiquidArgon
/mydet/isCalHomogeneous 0
#
/mydet/isUnitInLambda 1
/mydet/absorberTotalLength 9.96
/mydet/calorimeterRadius 4.98
/mydet/activeLayerNumber 60
/mydet/readoutLayerNumber 20
/mydet/activeLayerSize 8.5
/mydet/isRadiusUnitInLambda 1
/mydet/radiusBinSize 0.1
/mydet/radiusBinNumber 10
#
/mydet/update
#
/run/beamOn 5000
#
