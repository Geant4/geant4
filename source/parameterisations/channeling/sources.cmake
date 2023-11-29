# - G4channeling module build definition

# Define the Geant4 Module.
geant4_add_module(G4channeling
   PUBLIC_HEADERS
     G4BaierKatkov.hh
     G4ChannelingFastSimCrystalData.hh
     G4ChannelingFastSimInterpolation.hh
     G4ChannelingFastSimModel.hh
     G4VChannelingFastSimCrystalData.hh
   SOURCES
     G4BaierKatkov.cc
     G4ChannelingFastSimCrystalData.cc
     G4ChannelingFastSimInterpolation.cc
     G4ChannelingFastSimModel.cc
     G4VChannelingFastSimCrystalData.cc)

 geant4_module_link_libraries(G4channeling
   PUBLIC
     G4geometrymng
     G4globman
     G4heprandom
     G4materials
     G4parameterisation
     G4track
     G4partman
   PRIVATE
     G4bosons
     G4navigation)