#----------------------------------------------------------------
# 11-Aug-2004 A.R.   
#
# This is the prototype of the Geant4 command script for the
# statistical acceptance suite, based on a Calorimeter setup.
# The user has to make the following choices (in practice, 
# commenting/uncommeting few lines below):
#  1) Choice of the  * Particle Type * :
#       mu-, mu+, e-, e+, gamma, 
#       pi-, pi+, kaon-, kaon+, kaon0L, 
#       neutron, proton.
#  2) Choice of the  * Beam Energy * :
#       1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 80,
#       100, 120, 150, 180, 200, 250, 300, 1000    GeV.
#  3) Choice of the  * Calorimeter Type * :
#       I) Absorber and Active materials (and whether it is
#          an homogeneous calorimeter or a sampling one):
#            Fe-Sci, Cu-Sci, Cu-LAr, Pb-LAr, W-LAr  (Sampling) 
#            PbWO4  (Homogeneous)
#      II) Dimension and Segmentation:
#          - Is the unit with which to express the dimension of
#            the absorber in lambdas (interaction lengths), or
#            in [mm]?
#          - Size of the absorber in the given unit.
#          - Number of active layers (in the case of PbWO4, i.e.
#            of a homogeneous calorimeter, this number is relevant
#            for the longitudinal shower analysis).
#          - Size of the active layer, in [mm].
#          Then, for the transverse shower analysis:
#          - Is the unit with which to express the radius bin
#            size in lambdas (interaction lengths) of the absorber,
#            or in [mm]?	
#          - Size of the first bin for the radius (transverse 
#            distance from the beam axis): the other ones will
#            have increasing widths (e.g. the second bin has
#            a width equal to two times the first bin size,
#            the third one has a width equal to three times
#            the first bin size, etc.).
#          - Number of bins for the radius.
#     III) Update Geometry : always necessary, leave it On!
#
#----------------------------------------------------------------
#
/random/setSavingFlag 1
#
/run/verbose 1 
/event/verbose 0 
/tracking/verbose 0 
#
#=======================  PARTICLE TYPE  ====================
#
#/gun/particle mu-
#/gun/particle mu+
#/gun/particle e-
#/gun/particle e+
#/gun/particle gamma
#/gun/particle pi-
#/gun/particle pi+
#/gun/particle kaon-
#/gun/particle kaon+
#/gun/particle kaon0L
#/gun/particle neutron
/gun/particle proton
#
#=======================  BEAM ENERGY  ====================
#
#/gun/energy    1 GeV
#/gun/energy    2 GeV
#/gun/energy    3 GeV
#/gun/energy    4 GeV
#/gun/energy    5 GeV
#/gun/energy    6 GeV
#/gun/energy    7 GeV
#/gun/energy    8 GeV
#/gun/energy    9 GeV
#/gun/energy   10 GeV
/gun/energy   20 GeV
#/gun/energy   30 GeV
#/gun/energy   40 GeV
#/gun/energy   50 GeV
#/gun/energy   60 GeV
#/gun/energy   80 GeV
#/gun/energy  100 GeV
#/gun/energy  120 GeV
#/gun/energy  150 GeV
#/gun/energy  180 GeV
#/gun/energy  200 GeV
#/gun/energy  250 GeV
#/gun/energy  300 GeV
#/gun/energy 1000 GeV
#
#=======================  CALORIMETER TYPE  ====================
#
#=== I) ABSORBER and ACTIVE MATERIALS; and SAMPLING/HOMOGENEOUS TYPE ===
#
#--- Iron - Scintillator
/mydet/absorberMaterial Iron
/mydet/activeMaterial Scintillator
/mydet/isCalHomogeneous 0
#
#--- Copper - Scintillator
#/mydet/absorberMaterial Copper
#/mydet/activeMaterial Scintillator
#/mydet/isCalHomogeneous 0
#
#--- Copper - LiquidArgon
#/mydet/absorberMaterial Copper
#/mydet/activeMaterial LiquidArgon
#/mydet/isCalHomogeneous 0
#
#--- Lead - LiquidArgon
#/mydet/absorberMaterial Lead
#/mydet/activeMaterial LiquidArgon
#/mydet/isCalHomogeneous 0
#
#--- Tungsten - LiquidArgon
#/mydet/absorberMaterial Tungsten
#/mydet/activeMaterial LiquidArgon
#/mydet/isCalHomogeneous 0
#
#--- PbWO4 ---
#/mydet/absorberMaterial PbWO4
#/mydet/activeMaterial PbWO4
#/mydet/isCalHomogeneous 1
#
#=== II) DIMENSION and SEGMENTATION ===
#
/mydet/isUnitInLambda 1
/mydet/absorberTotalLength 10.0
/mydet/activeLayerNumber 20
/mydet/activeLayerSize 4.0
/mydet/isRadiusUnitInLambda 1
/mydet/radiusBinSize 0.1
/mydet/radiusBinNumber 10
#
#=== III) UPDATE GEOMETRY : leave it always ON ! ===
#
/mydet/update
#
#=======================  NUMBER OF EVENTS  ====================
#
#/run/beamOn    10
#/run/beamOn   100
#/run/beamOn  1000
/run/beamOn  5000
#/run/beamOn 10000
#

