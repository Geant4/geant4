#------------------------------------------------------------
# CMS electromagnetic crystal: homogeneous PbWO4.
#                              23 cm length.
#
# CMS hadronic barrel: sampling Cu-Sci.
#                      25 layers: 6 cm (Cu) - 4 mm (Sci).
#                      (total about 10 lambda).
#
# Radius of the cylinder: 75 cm
#                         (diameter about 10 lambda)
#------------------------------------------------------------
#
###/random/resetEngineFrom start.rndm
/random/setSavingFlag 1
#
/run/verbose 1 
/event/verbose 0 
/tracking/verbose 0 
#
#/gun/particle mu-
#/gun/particle e-
/gun/particle pi-
#/gun/particle proton
#
/gun/energy 5 GeV
#
/mydet/trackerLength 0.0
#
/mydet/emAbsorberMaterial PbWO4
/mydet/isEmCalHomogeneous 1
/mydet/emAbsorberTotalLength 230.0
#
/mydet/hadAbsorberMaterial Copper
/mydet/hadActiveMaterial Scintillator
/mydet/isHadCalHomogeneous 0
/mydet/hadAbsorberTotalLength 1500.0
/mydet/hadActiveLayerNumber 25
/mydet/hadActiveLayerSize 4.0
#
/mydet/muonLength 0.0
#
/mydet/detectorRadius 750.0
#
/mydet/update
#
/run/beamOn 100
#
