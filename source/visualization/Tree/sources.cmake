# - G4Tree module build definition

# Define the Geant4 Module.
geant4_add_module(G4Tree
  PUBLIC_HEADERS
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
    G4VTreeViewer.cc)

geant4_module_link_libraries(G4Tree
  PUBLIC
    G4modeling
    G4intercoms
    G4vis_management
  PRIVATE
    G4graphics_reps
    G4geometrymng
    G4materials
    G4navigation
    G4globman
    G4detector)
