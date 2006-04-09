#----------------------------------------------------------------
# Last update: 09-Apr-2006
#
# There are three proposals for the ILC electromagnetic 
# calorimeter, and five for the ILC hadronic calorimeter.
# We choose, to simplify, the proposal for ECAL and the
# one for HCAL which coincides, a part on the number of
# layers: 
#                Pb (4.0 mm) - Sci (1 mm)
#
# For  Pb :         X0          Lambda     Lambda/X0
#                  0.5612 cm    17.092 cm    30.46
# so:
#                  4.0 mm = 0.713 X0 = 0.0234 lambda
#
# For the ECAL part, the proposal consists of   38  layers;
# for the HCAL part, the proposal consists of  130  layers.
#
# The use the following configuration:  168  active layers ,
# each layer made of precisely  4.0 mm Pb  and  1.0 mm Sci .
# We use  24 readout layers : so the ECAL part would be
# the first 5 readout layers, plus half of the 6th; the 
# rest is the hadronic part.
# For the lateral profile, we use  30 bins  (so more than usual
# because we want to study electromagnetic and hadronic shower
# shapes at the same time), with  X0(Pb)/4 = 1.403 mm  as size 
# of the first bin.
#
#----------------------------------------------------------------
#
/random/setSavingFlag 1
#
/run/verbose 1 
/event/verbose 0 
/tracking/verbose 0 
#
/gun/particle e-
#/gun/particle pi-
#
/gun/energy 100 GeV
#
/mydet/absorberMaterial Lead
/mydet/activeMaterial Scintillator
/mydet/isCalHomogeneous 0
#
/mydet/isUnitInLambda 0
/mydet/absorberTotalLength 672.0
/mydet/calorimeterRadius 336.0
/mydet/activeLayerNumber 168
/mydet/readoutLayerNumber 24
/mydet/activeLayerSize 1.0
/mydet/isRadiusUnitInLambda 0
/mydet/radiusBinSize 1.403
/mydet/radiusBinNumber 30
#
/mydet/update
#
/run/beamOn 5000
#

