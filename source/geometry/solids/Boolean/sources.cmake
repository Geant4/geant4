#------------------------------------------------------------------------------
# Module : G4geomBoolean
# Package: Geant4.src.G4geometry.G4geomBoolean
#------------------------------------------------------------------------------

#
# Define the Geant4 Module.
#
geant4_define_module(NAME G4geomBoolean
  HEADERS
    G4BooleanSolid.hh
    G4BooleanSolid.icc
    G4DisplacedSolid.hh
    G4IntersectionSolid.hh
    G4MultiUnion.hh
    G4ScaledSolid.hh
    G4SubtractionSolid.hh
    G4UnionSolid.hh
  SOURCES
    G4BooleanSolid.cc
    G4DisplacedSolid.cc
    G4IntersectionSolid.cc
    G4MultiUnion.cc
    G4ScaledSolid.cc
    G4SubtractionSolid.cc
    G4UnionSolid.cc
  GRANULAR_DEPENDENCIES
    G4geometrymng
    G4globman
    G4graphics_reps
    G4intercoms
    G4volumes
    G4csg
    G4specsolids
  GLOBAL_DEPENDENCIES
    G4global
    G4graphics_reps
    G4intercoms
  LINK_LIBRARIES
    ${VECGEOM_LIBRARIES}
)

# List any source specific properties here
