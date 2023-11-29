# - G4interfaces category build
include(core/sources.cmake)
include(implementation/sources.cmake)

geant4_add_category(G4interfaces MODULES G4UIcore G4UIimplementation)
