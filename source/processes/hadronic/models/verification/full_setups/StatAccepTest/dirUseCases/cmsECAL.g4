#----------------------------------------------------------------
# Last update: 24-Jan-2006
#
# CMS ECAL : made of  2 cm x 2 cm x 23 cm  PbWO4 crystals.
#  
# For PbWO4 :        X0          Lambda     Lambda/X0
#                   0.89 cm      22.4 cm      25.17
# so:
#                    23 cm = 25.8 X0 = 1.03 lambda
#
# The configuration which we used is a cube of  23 cm  length 
# with only 1 layer: notice that we cannot avoid to use an
# active layer, even if we do not need it. However, we can
# get the energy in the "absorber" layer, which is what we 
# are interested in, as:
#                         EDEP_CAL - EDEP_ACT
#
# where the two ntuple's variables EDEP_CAL and EDEP_ACT
# are, respectively, the total energy deposit in the whole
# calorimeter (passive + active), and the energy deposited
# in the active layer.
# It is not possible to look at the transverse profile, 
# because with one layer it gives energy in the active layer
# only, not in the "absorber" one in which we are interested.
# Notice that, being not interested to the "fictitious" 
# active layer, we do not specify anything about it, taking
# whatever default is provided. The same for the transverse
# profile.
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
/mydet/absorberMaterial PbWO4
#
/mydet/isUnitInLambda 0
/mydet/absorberTotalLength 230
/mydet/activeLayerNumber 1
/mydet/readoutLayerNumber 1
#
###/mydet/setField 4.0 tesla
#
/mydet/update
#
/run/beamOn 5000
#