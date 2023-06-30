# - G4geomBoolean module build definition

# Define the Geant4 Module.
geant4_add_module(G4geomBoolean
  PUBLIC_HEADERS
    G4BooleanSolid.hh
    G4BooleanSolid.icc
    G4DisplacedSolid.hh
    G4IntersectionSolid.hh
    G4MultiUnion.hh
    G4ScaledSolid.hh
    G4SubtractionSolid.hh
    G4UnionSolid.hh
    G4VBooleanProcessor.hh
  SOURCES
    G4BooleanSolid.cc
    G4DisplacedSolid.cc
    G4IntersectionSolid.cc
    G4MultiUnion.cc
    G4ScaledSolid.cc
    G4SubtractionSolid.cc
    G4UnionSolid.cc)

geant4_module_link_libraries(G4geomBoolean
  PUBLIC G4geometrymng G4hepgeometry G4globman G4specsolids ${VECGEOM_LIBRARIES}
  PRIVATE G4graphics_reps G4heprandom)
