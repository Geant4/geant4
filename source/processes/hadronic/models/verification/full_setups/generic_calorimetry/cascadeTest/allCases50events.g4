#----------------------------------------------------------------
# All 735 = 7 calorimeter type x 7 particle types x 15 energies 
# cases, with 50 events each --> 36,750 events.
# This corrsponds to 1/100 of the wanted 5,000 events for each
# case.
#----------------------------------------------------------------
#
/random/setSavingFlag 1
#
/run/verbose 0 
/event/verbose 0 
/tracking/verbose 0 
#
#--- Iron - Scintillator
/mydet/absorberMaterial Iron
/mydet/activeMaterial Scintillator
#
#    --- pi- ---
/gun/particle pi-
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#    --- pi+ ---
/gun/particle pi+
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#    --- kaon- ---
/gun/particle kaon-
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#    --- kaon+ ---
/gun/particle kaon+
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#    --- kaon0L ---
/gun/particle kaon0L
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#    --- proton ---
/gun/particle proton
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#    --- neutron ---
/gun/particle neutron
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#--- Copper - Scintillator
/mydet/absorberMaterial Copper
/mydet/activeMaterial Scintillator
#
#    --- pi- ---
/gun/particle pi-
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#    --- pi+ ---
/gun/particle pi+
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#    --- kaon- ---
/gun/particle kaon-
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#    --- kaon+ ---
/gun/particle kaon+
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#    --- kaon0L ---
/gun/particle kaon0L
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#    --- proton ---
/gun/particle proton
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#    --- neutron ---
/gun/particle neutron
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#--- Copper - LiquidArgon
/mydet/absorberMaterial Copper
/mydet/activeMaterial LiquidArgon
#
#    --- pi- ---
/gun/particle pi-
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#    --- pi+ ---
/gun/particle pi+
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#    --- kaon- ---
/gun/particle kaon-
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#    --- kaon+ ---
/gun/particle kaon+
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#    --- kaon0L ---
/gun/particle kaon0L
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#    --- proton ---
/gun/particle proton
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#    --- neutron ---
/gun/particle neutron
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#--- Lead - LiquidArgon
/mydet/absorberMaterial Lead
/mydet/activeMaterial LiquidArgon
#
#    --- pi- ---
/gun/particle pi-
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#    --- pi+ ---
/gun/particle pi+
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#    --- kaon- ---
/gun/particle kaon-
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#    --- kaon+ ---
/gun/particle kaon+
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#    --- kaon0L ---
/gun/particle kaon0L
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#    --- proton ---
/gun/particle proton
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#    --- neutron ---
/gun/particle neutron
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#--- PbWO4 ---
/mydet/absorberMaterial PbWO4
/mydet/activeMaterial PbWO4
#
#    --- pi- ---
/gun/particle pi-
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#    --- pi+ ---
/gun/particle pi+
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#    --- kaon- ---
/gun/particle kaon-
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#    --- kaon+ ---
/gun/particle kaon+
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#    --- kaon0L ---
/gun/particle kaon0L
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#    --- proton ---
/gun/particle proton
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#    --- neutron ---
/gun/particle neutron
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#--- Tungsten - Silicon
/mydet/absorberMaterial Tungsten
/mydet/activeMaterial Silicon
#
#    --- pi- ---
/gun/particle pi-
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#    --- pi+ ---
/gun/particle pi+
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#    --- kaon- ---
/gun/particle kaon-
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#    --- kaon+ ---
/gun/particle kaon+
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#    --- kaon0L ---
/gun/particle kaon0L
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#    --- proton ---
/gun/particle proton
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#    --- neutron ---
/gun/particle neutron
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#--- Tungsten - LiquidArgon
/mydet/absorberMaterial Tungsten
/mydet/activeMaterial LiquidArgon
#
#    --- pi- ---
/gun/particle pi-
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#    --- pi+ ---
/gun/particle pi+
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#    --- kaon- ---
/gun/particle kaon-
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#    --- kaon+ ---
/gun/particle kaon+
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#    --- kaon0L ---
/gun/particle kaon0L
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#    --- proton ---
/gun/particle proton
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
#    --- neutron ---
/gun/particle neutron
#
/gun/energy 3.6 GeV
/run/beamOn 50
#
/gun/energy 5 GeV
/run/beamOn 50
#
/gun/energy 10 GeV
/run/beamOn 50
#
/gun/energy 20 GeV
/run/beamOn 50
#
/gun/energy 30 GeV
/run/beamOn 50
#
/gun/energy 40 GeV
/run/beamOn 50
#
/gun/energy 50 GeV
/run/beamOn 50
#
/gun/energy 60 GeV
/run/beamOn 50
#
/gun/energy 80 GeV
/run/beamOn 50
#
/gun/energy 100 GeV
/run/beamOn 50
#
/gun/energy 120 GeV
/run/beamOn 50
#
/gun/energy 150 GeV
/run/beamOn 50
#
/gun/energy 180 GeV
/run/beamOn 50
#
/gun/energy 200 GeV
/run/beamOn 50
#
/gun/energy 300 GeV
/run/beamOn 50
#
##--- Uranium - LiquidArgon
###/mydet/absorberMaterial Uranium
##/mydet/activeMaterial LiquidArgon
##



