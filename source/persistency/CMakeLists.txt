# - G4persistency category build

set(G4persistency_COMPONENTS
  ascii/sources.cmake
  mctruth/sources.cmake)

if(GEANT4_USE_GDML)
  list(APPEND G4persistency_COMPONENTS gdml/sources.cmake)
endif()

geant4_global_library_target(NAME G4persistency COMPONENTS ${G4persistency_COMPONENTS})
