#----------------------------------------------------------------
# Last update: 10-Jul-2007
#
# ATLAS Barrel LAr electromagnetic calorimeter:
#                Pb (1.53 mm) - LAr (4.2 mm)
#
# For  Pb :         X0          Lambda     Lambda/X0
#                  0.5612 cm    17.092 cm    30.46
# so:
#                  1.53 mm = 0.27 X0 = 0.0090 lambda
#
# For  LAr :        X0          Lambda     Lambda/X0
#                  14.0036 cm   85.7066 cm    6.12
# so:
#                  4.2 mm = 0.03 X0 = 0.0049 lambda
#
# Hence a whole layer (passive+active): 0.30 X0 = 0.0139 lambda
#
# To get as close as possible to the total thickness of the
# real calorimeter, 27.6 X0, we use the following configuration:  
#  92  active layers ,
# each layer made of precisely  1.53 mm Pb  and  4.2 mm LAr .
# We use  23 readout layers  (between the divisors of the
# number of active layers (92), 23 is the closest to the
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
# The real setup is made of 3 longitudinal segmentation of 
# thicknesses: 4.7 X0 (front); 18.1 X0 (middle); 4.8 X0 (back).
# A reasonable way to emulate these three segmentations is the
# following: 
#   o  front :  first  16 active layers, so first   4 readout ones;
#   o  middle : middle 60 active layers, so middle 15 readout ones;
#   o  back :   last   16 active layers, so last    4 readout ones.
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
/mydet/activeMaterial LiquidArgon
/mydet/isCalHomogeneous 0
#
/mydet/isUnitInLambda 0
/mydet/absorberTotalLength 140.76
/mydet/calorimeterRadius 854.6
/mydet/activeLayerNumber 92
/mydet/readoutLayerNumber 23
/mydet/activeLayerSize 4.2
/mydet/isRadiusUnitInLambda 0
/mydet/radiusBinSize 1.403
/mydet/radiusBinNumber 10
#
###/mydet/setField 2.0 tesla
#
/mydet/update
#
/run/beamOn 5000
#

