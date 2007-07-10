#----------------------------------------------------------------
# Last update: 10-Jul-2007
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
#
# As radius of the cylindrical calorimeter, we take a large
# one, corresponding to the usual  5 lambda radius  
# ("lambda" of the absorber material, Pb) used for
# hadronic calorimeter: in practice, for electromagnetic
# showers, this corresponds to an infinite wide calorimeter,
# without lateral leakage. 
#
# For the binning of lateral profile, we use (as usual) 
# 10 bins, but with X0(Pb)/4 = 1.403 mm as size of the first
# bin, in such a way to have the proper "granularity" to
# study an electromagnetic shower. This implies that the
# last lateral bin is very big: for electromagnetic showers
# it should anyhow being almost empty, but for hadronic 
# showers it should get a large fraction of visible energy.
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
/mydet/calorimeterRadius 854.6
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

