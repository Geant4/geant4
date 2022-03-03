# - G4VRML module build definition

# Define the Geant4 Module.
geant4_add_module(G4VRML
  PUBLIC_HEADERS
    G4VRML2File.hh
    G4VRML2FileSceneHandler.hh
    G4VRML2FileViewer.hh
  SOURCES
    G4VRML2File.cc
    G4VRML2FileSceneHandler.cc
    G4VRML2FileViewer.cc
    G4VRML2SceneHandlerFunc.icc)

geant4_module_link_libraries(G4VRML
  PUBLIC
    G4geometrymng
    G4vis_management
    G4globman
  PRIVATE
    G4csg
    G4graphics_reps
    G4modeling)
