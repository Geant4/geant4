# - G4FR module build definition

# Define the Geant4 Module.
geant4_add_module(G4FR
  PUBLIC_HEADERS
    G4DAWNFILE.hh
  PRIVATE_HEADERS
    G4DAWNFILESceneHandler.hh
    G4DAWNFILEViewer.hh
    G4FRConst.hh
    G4FRSceneFunc.icc
    G4FRofstream.hh
    G4VisFeaturesOfDAWNFILE.hh
  SOURCES
    G4DAWNFILE.cc
    G4DAWNFILESceneHandler.cc
    G4DAWNFILEViewer.cc
    G4VisFeaturesOfDAWNFILE.cc)

geant4_module_link_libraries(G4FR
  PUBLIC
    G4vis_management
  PRIVATE
    G4csg
    G4geometrymng
    G4globman
    G4graphics_reps
    G4hepgeometry
    G4modeling)

