#--------------------------------------------------------------------
# Last update: 16-Nov-2007
#
# ATLAS TileCal at 90 degrees:  Fe(14mm)-Sci(3mm) .
# 
# We use 120 layers, which corresponds to about 10 lambda of
# absorber (Fe), and 20 readout layers.
#
# For  Fe :         X0          Lambda     Lambda/X0
#                  1.7585 cm    16.760 cm     9.53    
#  so:
#       1 Iron layer = 1.4 cm = 0.0835 lambda = 0.796 X0
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
/mydet/isUnitInLambda 0
/mydet/absorberTotalLength 1680.0
/mydet/calorimeterRadius 840.0
/mydet/activeLayerNumber 120
/mydet/readoutLayerNumber 20
/mydet/activeLayerSize 3.0
/mydet/isRadiusUnitInLambda 1
/mydet/radiusBinSize 0.1
/mydet/radiusBinNumber 10
#
/mydet/update
#
/run/beamOn 5000
#
