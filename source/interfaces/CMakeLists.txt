# - G4interfaces category build
include(basic/sources.cmake)
include(common/sources.cmake)

geant4_add_category(G4interfaces MODULES G4UIbasic G4UIcommon)
