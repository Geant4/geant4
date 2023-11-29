# - G4parmodels category build
include(channeling/sources.cmake)
include(gflash/sources.cmake)
geant4_add_category(G4parmodels MODULES G4channeling G4gflash)
