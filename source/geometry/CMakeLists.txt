# - G4geometry category builds
add_definitions(-DG4GEOM_ALLOC_EXPORT)

include(biasing/sources.cmake)
include(divisions/sources.cmake)
include(magneticfield/sources.cmake)
include(management/sources.cmake)
include(navigation/sources.cmake)
include(solids/Boolean/sources.cmake)
include(solids/CSG/sources.cmake)
include(solids/specific/sources.cmake)
include(volumes/sources.cmake)

geant4_add_category(G4geometry MODULES
  G4geombias
  G4geomdivision
  G4magneticfield
  G4geometrymng
  G4navigation
  G4geomBoolean
  G4csg
  G4specsolids
  G4volumes)
