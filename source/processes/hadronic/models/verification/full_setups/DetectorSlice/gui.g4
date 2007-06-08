#
# This file permits to customize, with commands,
# the menu bar of the G4UIXm sessions.
# It has no effect with G4UIterminal.
#
# --- File --- 
/gui/addMenu file File
/gui/addButton file Continue continue
/gui/addButton file Exit "exit"
#
# --- B-field along z ---
/gui/addMenu Bfield BField
/gui/addButton Bfield 0-Tesla  "/mydet/setField 0 tesla" 
/gui/addButton Bfield 1-Tesla  "/mydet/setField 1 tesla" 
/gui/addButton Bfield 2-Tesla  "/mydet/setField 2 tesla" 
/gui/addButton Bfield 3-Tesla  "/mydet/setField 3 tesla" 
/gui/addButton Bfield 4-Tesla  "/mydet/setField 4 tesla"
# 
# --- Materials ---
/gui/addMenu absorber Absorber
/gui/addButton absorber Iron     "/mydet/absorberMaterial Iron"
/gui/addButton absorber Copper   "/mydet/absorberMaterial Copper"
/gui/addButton absorber Tungsten "/mydet/absorberMaterial Tungsten"
/gui/addButton absorber PbWO4    "/mydet/absorberMaterial PbWO4"
/gui/addButton absorber Lead     "/mydet/absorberMaterial Lead"
/gui/addButton absorber Uranium  "/mydet/absorberMaterial Uranium"
#
/gui/addMenu active Active
/gui/addButton active Scintillator "/mydet/activeMaterial Scintillator"
/gui/addButton active LiquidArgon  "/mydet/activeMaterial LiquidArgon"
/gui/addButton active PbWO4        "/mydet/activeMaterial PbWO4"
/gui/addButton active Silicon      "/mydet/activeMaterial Silicon"
/gui/addButton active Quartz       "/mydet/activeMaterial Quartz"
#
# --- Particle ---
/gui/addMenu particle Particle
/gui/addButton particle geantino         "/gun/particle geantino"
/gui/addButton particle chargedgeantino  "/gun/particle chargedgeantino"
/gui/addButton particle muon-            "/gun/particle mu-"
/gui/addButton particle electron-        "/gun/particle e-"
/gui/addButton particle pion-            "/gun/particle pi-"
/gui/addButton particle pion+            "/gun/particle pi+"
/gui/addButton particle kaon-            "/gun/particle kaon-"
/gui/addButton particle kaon+            "/gun/particle kaon+"
/gui/addButton particle kaon0L           "/gun/particle kaon0L"
/gui/addButton particle proton           "/gun/particle proton"
/gui/addButton particle neutron          "/gun/particle neutron"
#
# --- Energy ---
/gui/addMenu energy Energy
/gui/addButton energy   3.6-GeV  "/gun/energy  3.6 GeV"
/gui/addButton energy    5-GeV  "/gun/energy   5 GeV"
/gui/addButton energy   10-GeV  "/gun/energy  10 GeV"
/gui/addButton energy   20-GeV  "/gun/energy  20 GeV"
/gui/addButton energy   30-GeV  "/gun/energy  30 GeV"
/gui/addButton energy   40-GeV  "/gun/energy  40 GeV"
/gui/addButton energy   50-GeV  "/gun/energy  50 GeV"
/gui/addButton energy   60-GeV  "/gun/energy  60 GeV"
/gui/addButton energy   80-GeV  "/gun/energy  80 GeV"
/gui/addButton energy  100-GeV  "/gun/energy 100 GeV"
/gui/addButton energy  120-GeV  "/gun/energy 120 GeV"
/gui/addButton energy  150-GeV  "/gun/energy 150 GeV"
/gui/addButton energy  180-GeV  "/gun/energy 180 GeV"
/gui/addButton energy  200-GeV  "/gun/energy 200 GeV"
/gui/addButton energy  300-GeV  "/gun/energy 300 GeV"
#
# --- NumEvents ---
/gui/addMenu numEvents NumEvents
/gui/addButton numEvents    1   "/run/beamOn    1"
/gui/addButton numEvents    2   "/run/beamOn    2"
/gui/addButton numEvents    5   "/run/beamOn    5"
/gui/addButton numEvents   10   "/run/beamOn   10"
/gui/addButton numEvents   20   "/run/beamOn   20"
/gui/addButton numEvents   50   "/run/beamOn   50"
/gui/addButton numEvents  100   "/run/beamOn  100"
/gui/addButton numEvents 5000   "/run/beamOn  5000"
#
