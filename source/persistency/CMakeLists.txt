# - G4persistency category build
include(ascii/sources.cmake)
include(mctruth/sources.cmake)

geant4_add_category(G4persistency MODULES G4geomtext G4mctruth)

if(GEANT4_USE_GDML)
  include(gdml/sources.cmake)
  geant4_category_modules(G4persistency G4gdml)
endif()
