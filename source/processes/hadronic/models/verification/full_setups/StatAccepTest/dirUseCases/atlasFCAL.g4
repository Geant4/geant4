#----------------------------------------------------------------
# Last update: 14-Jun-2006
#
# This is a very simplfied version of the ATLAS Forward calorimeter,
# (in particular the second hadronic module)
# which is made of 5 mm diameter tubes (with axis along the beam
# direction), containing Tungsten rods of  4.5 mm diameter , 
# with  0.5 mm LAr  gaps. 
#
# For W : X0 = 0.35 cm  and  Lambda = 9.5855 cm
#
# so      4.5 mm = 1.29 X0 = 0.046946 lambda
#  
# The closest configuration with an overall thickness of about
# the usual 10 lambda (more precisely, 9.85864 lambda) is the 
# following:  210  active layers , 
# which corresponds to  4.5 mm , i.e. 0.046946 lambda  
# per passive layer.           
# We use  21  readout layers.
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
/gun/energy  100 GeV
#
/mydet/absorberMaterial Tungsten
/mydet/activeMaterial LiquidArgon
/mydet/isCalHomogeneous 0
#
/mydet/isUnitInLambda 1
/mydet/absorberTotalLength 9.85864
/mydet/calorimeterRadius 5.0
/mydet/activeLayerNumber 210
/mydet/readoutLayerNumber 21
/mydet/activeLayerSize 0.5
/mydet/isRadiusUnitInLambda 1
/mydet/radiusBinSize 0.1
/mydet/radiusBinNumber 10
#
###/mydet/setField 4.0 tesla
#
/mydet/update
#
/run/beamOn  5000
#

