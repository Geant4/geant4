#----------------------------------------------------------------
# Last update: 09-Apr-2006
#
# LHCb electromagnetic calorimeter:  Pb (2.0 mm) - Sci (4 mm)
#
# For  Pb :         X0          Lambda     Lambda/X0
#                  0.5612 cm    17.092 cm    30.46
# so:
#                  2.0 mm = 0.356 X0 = 0.0117 lambda
#
# We use the following configuration:  70  active layers ,
# each layer made of precisely  2.0 mm Pb  and  4.0 mm Sci .
# We use  14 readout layers  (between the divisors of the
# number of active layers (70), 14 is the closest to the
# usual (20) number of readout layers we are used to consider
# for the longitudinal profile).
# For the lateral profile, we use 10 bins (as usual), with
# X0(Pb)/4 = 1.403 mm as size of the first bin.
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
/mydet/absorberTotalLength 140.0
/mydet/calorimeterRadius 70.0
/mydet/activeLayerNumber 70
/mydet/readoutLayerNumber 14
/mydet/activeLayerSize 4.0
/mydet/isRadiusUnitInLambda 0
/mydet/radiusBinSize 1.403
/mydet/radiusBinNumber 10
#
/mydet/update
#
/run/beamOn 5000
#

