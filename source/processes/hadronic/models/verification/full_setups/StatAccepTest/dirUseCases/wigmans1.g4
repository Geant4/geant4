#----------------------------------------------------------------
# Last update: 09-Apr-2006
#
# First of the two Wigmans' proposed challenges for Geant4, 
# based on the paper: NIM A262 (1987) 229.
#
#               Pb (10.0 mm) - Sci (2.5 mm)
#
# For  Pb :         X0          Lambda     Lambda/X0
#                  0.5612 cm    17.092 cm    30.46
# so:
#                  10.0 mm = 1.78 X0 = 0.0585 lambda
#
# Although the actual calorimeter was 1 (EM) + 4 (HAD) lambda thick,
# we will use as usual (roughly) 10 lambda, to avoid problems with
# leakages. To have about 10 lambdas, we need roughly 171 layers.
# Luckly, this is divisible by 9, so we have 19 readout layers,
# although the shower profile is not treated in the paper.
# NB) To get the EM part, corresponding to 1 Lambda, it is enough
#     roughly to sum the first two readout layers.
#
# Beam energies: 3 , 5 , 7 , 10 , 20 , 30 , 50 , 75  GeV .
#
# Energy resolution: 
#   - Electrons:  sigma_e / <E_e> = A/sqrt(E)  "quadratic sum"  B
#                 A = (23.5 +/- 0.2)%  and  B = (1.2 +/- 0.2)%
#   - Pion-:      sigma_h / < E_h > = (44.2 +/- 1.3)% / sqrt(E)
#
# Response ratio:
#   - e/h : almost energy independent for beam energies above
#           10 GeV :    e/h = 1.05 +/- 0.04   for E > 10 GeV .
#   - e/mip (using muons of energies 5, 10, 50 GeV in the
#            electromagnetic section, i.e. first lambda of the
#            calorimeter):  e/mip = 0.67 +/- 0.032 .
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
###/gun/particle mu-
###/gun/particle pi-
#
###/gun/energy 3 GeV
###/gun/energy 5 GeV
###/gun/energy 7 GeV
/gun/energy 10 GeV
###/gun/energy 20 GeV
###/gun/energy 30 GeV
###/gun/energy 50 GeV
###/gun/energy 75 GeV
#
/mydet/absorberMaterial Lead
/mydet/activeMaterial Scintillator
/mydet/isCalHomogeneous 0
#
/mydet/isUnitInLambda 1
/mydet/absorberTotalLength 10.0035
/mydet/calorimeterRadius 5.00175
/mydet/activeLayerNumber 171
/mydet/readoutLayerNumber 19
/mydet/activeLayerSize 2.5
/mydet/isRadiusUnitInLambda 1
/mydet/radiusBinSize 0.1
/mydet/radiusBinNumber 10
#
/mydet/update
#
/run/beamOn 5000
#

