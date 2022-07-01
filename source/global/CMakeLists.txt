# - G4global category build(s)
add_definitions(-DG4GLOB_ALLOC_EXPORT)

include(HEPGeometry/sources.cmake)
include(HEPNumerics/sources.cmake)
include(HEPRandom/sources.cmake)
include(management/sources.cmake)

geant4_add_category(G4global MODULES 
  G4hepgeometry
  G4hepnumerics
  G4heprandom
  G4globman )