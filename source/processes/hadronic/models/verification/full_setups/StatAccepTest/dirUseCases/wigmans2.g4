#----------------------------------------------------------------
# Last update: 09-Apr-2006
#
# Second of the two Wigmans' proposed challenges for Geant4, 
# based on the paper: NIM A274 (1989) 134.
#
#               Pb (5.0 mm) - Sci (5.0 mm)
#
# For  Pb :         X0          Lambda     Lambda/X0
#                  0.5612 cm    17.092 cm    30.46
# so:
#                  5.0 mm = 0.89 X0 = 0.02925 lambda
#
# Although the actual calorimeter was 4,2 lambda thick,
# we will use as usual (roughly) 10 lambda, to avoid problems with
# leakages. To have about 10 lambdas, we need roughly 342 layers.
# Luckly, this is divisible by 18, so we have 19 readout layers,
# although the shower profile is not treated in the paper.
#
# Beam energies: 3 , 5 , 7 , 8.75  GeV .
# However, we consider here the same beam energies as in the first
# of Wigmans' test cases: 3 , 5 , 7 , 10 , 20 , 30 , 50 , 75  GeV .
#
# Energy resolution: 
#   - Electrons:  sigma_e / <E_e> = A/sqrt(E)  "quadratic sum"  B
#                 A = (12.7 +/- ???)%  and  B = (1.2 +/- ???)%
#                 (Wigmans reported: sigma_e / <E_e> = 14% / sqrt(E) )
#   - Pion-:      not reported in the paper, nor explicitly by
#                 Wigmans, but:
#                     sigma_h / <E_h> = (>>44%) / sqrt(E)
#
# Response ratio:
#   - e/h :  @3 GeV      e/h = 1.36 +/- 0.06
#            @5 GeV      e/h = 1.33 +/- 0.06
#            @6 GeV      e/h = 1.35 +/- 0.06
#            @8.75 GeV   e/h = 1.34 +/- 0.06
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
/mydet/activeLayerNumber 342
/mydet/readoutLayerNumber 19
/mydet/activeLayerSize 5.0
/mydet/isRadiusUnitInLambda 1
/mydet/radiusBinSize 0.1
/mydet/radiusBinNumber 10
#
/mydet/update
#
/run/beamOn 5000
#

