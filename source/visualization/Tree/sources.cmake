#------------------------------------------------------------------------------
# Module : G4Tree
# Package: Geant4.src.G4visualization.G4Tree
#------------------------------------------------------------------------------

#
# Define the Geant4 Module.
#
geant4_define_module(NAME G4Tree
  HEADERS
    G4ASCIITree.hh
    G4ASCIITreeMessenger.hh
    G4ASCIITreeSceneHandler.hh
    G4ASCIITreeViewer.hh
    G4VTree.hh
    G4VTreeSceneHandler.hh
    G4VTreeSceneHandler.icc
    G4VTreeViewer.hh
  SOURCES
    G4ASCIITree.cc
    G4ASCIITreeMessenger.cc
    G4ASCIITreeSceneHandler.cc
    G4ASCIITreeViewer.cc
    G4VTree.cc
    G4VTreeSceneHandler.cc
    G4VTreeViewer.cc
  GRANULAR_DEPENDENCIES
    G4csg
    G4detector
    G4geometrymng
    G4globman
    G4graphics_reps
    G4hits
    G4intercoms
    G4materials
    G4modeling
    G4navigation
    G4partman
    G4specsolids
    G4track
    G4vis_management
    G4volumes
  GLOBAL_DEPENDENCIES
    G4digits_hits
    G4geometry
    G4global
    G4graphics_reps
    G4intercoms
    G4materials
    G4modeling
    G4particles
    G4track
    G4vis_management
)

# List any source specific properties here
