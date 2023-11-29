# - G4persistency category build
include(ascii/sources.cmake)
geant4_add_category(G4geomtext MODULES G4geomtext)

include(mctruth/sources.cmake)
geant4_add_category(G4mctruth MODULES G4mctruth)

if(GEANT4_USE_GDML)
  include(gdml/sources.cmake)
  geant4_add_category(G4gdml MODULES G4gdml)
endif()
