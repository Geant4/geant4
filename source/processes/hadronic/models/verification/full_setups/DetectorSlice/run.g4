#
###/random/resetEngineFrom start.rndm
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
/gun/particle pi-
#/gun/particle pi+
#/gun/particle kaon-
#/gun/particle kaon+
#/gun/particle kaon0L
#/gun/particle neutron
#/gun/particle proton
#/gun/particle anti_proton
#/gun/particle anti_neutron
#/gun/particle deuteron 
#/gun/particle triton 
#/gun/particle alpha 
#
#=======================  BEAM ENERGY  ====================
#
/gun/energy  10 GeV
#
#=======================  CALORIMETER TYPE  ====================
#
# --- Tracker ---
/mydet/trackerMaterial Silicon
###/mydet/trackerMaterial LiquidArgon
###/mydet/trackerMaterial Scintillator
#
# --- EM Calorimeter ---
/mydet/emAbsorberMaterial Lead
###/mydet/emAbsorberMaterial PbWO4
###/mydet/emAbsorberMaterial Iron
###/mydet/emAbsorberMaterial Copper
###/mydet/emAbsorberMaterial Tungsten
###/mydet/emAbsorberMaterial Uranium
#
/mydet/emActiveMaterial LiquidArgon
###/mydet/emActiveMaterial Scintillator
###/mydet/emActiveMaterial PbWO4
#
/mydet/isEmCalHomogeneous 0
###/mydet/isEmCalHomogeneous 1
#
# --- HAD Calorimeter ---
/mydet/hadAbsorberMaterial Iron
###/mydet/hadAbsorberMaterial Copper
###/mydet/hadAbsorberMaterial Tungsten
###/mydet/hadAbsorberMaterial Uranium
###/mydet/hadAbsorberMaterial Lead
###/mydet/hadAbsorberMaterial PbWO4
#
###/mydet/hadActiveMaterial LiquidArgon
/mydet/hadActiveMaterial Scintillator
###/mydet/hadActiveMaterial PbWO4
#
/mydet/isHadCalHomogeneous 0
###/mydet/isHadCalHomogeneous 1
#
# --- Muon detector ---
/mydet/muonMaterial Iron
###/mydet/muonMaterial Copper
###/mydet/muonMaterial Tungsten
###/mydet/muonMaterial Lead
###/mydet/muonMaterial Uranium
#
#=== II) DIMENSION and SEGMENTATION ===
#
/mydet/trackerLength 10.0
#
/mydet/emAbsorberTotalLength 200.0
/mydet/emActiveLayerNumber 50
/mydet/emActiveLayerSize 1.0
#
/mydet/hadAbsorberTotalLength 2000.0
/mydet/hadActiveLayerNumber 50
/mydet/hadActiveLayerSize 4.0
#
/mydet/muonLength 1000.0
#
###/mydet/detectorRadius 1000.0
/mydet/detectorRadius 2000.0
#
#=== III) UPDATE GEOMETRY : leave it always ON ! ===
#
/mydet/update
#
#=======================  NUMBER OF EVENTS  ====================
#
#/run/beamOn    10
/run/beamOn   100
#/run/beamOn  1000
#/run/beamOn  5000
#/run/beamOn 10000
#

