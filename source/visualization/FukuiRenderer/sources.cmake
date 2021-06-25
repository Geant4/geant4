# - G4FR module build definition

# Module has optional sources
# List those always built
set(G4VIS_DAWN_MODULE_HEADERS
  G4DAWNFILE.hh
  G4DAWNFILESceneHandler.hh
  G4DAWNFILEViewer.hh
  G4FRConst.hh
  G4FRFeatures.hh
  G4FRSceneFunc.icc
  G4FRofstream.hh
  G4VisFeaturesOfDAWNFILE.hh
  G4VisFeaturesOfFukuiRenderer.hh)

set(G4VIS_DAWN_MODULE_SOURCES
  G4DAWNFILE.cc
  G4DAWNFILESceneHandler.cc
  G4DAWNFILEViewer.cc
  G4VisFeaturesOfDAWNFILE.cc
  G4VisFeaturesOfFukuiRenderer.cc)

# DAWN Network driver only built if user selected
if(GEANT4_USE_NETWORKDAWN)
  list(APPEND G4VIS_DAWN_MODULE_HEADERS
    G4FRClientServer.hh
    G4FukuiRenderer.hh
    G4FukuiRendererSceneHandler.hh
    G4FukuiRendererViewer.hh)

  list(APPEND G4VIS_DAWN_MODULE_SOURCES
    G4FRClientServer.cc
    G4FukuiRenderer.cc
    G4FukuiRendererSceneHandler.cc
    G4FukuiRendererViewer.cc)

  # To activate the Fukui Renderer Network component, we need an
  # extra compile definition
  add_definitions(-DG4VIS_BUILD_DAWN_DRIVER)
endif()

# Define the Geant4 Module.
geant4_add_module(G4FR
  PUBLIC_HEADERS
    ${G4VIS_DAWN_MODULE_HEADERS}
  SOURCES
    ${G4VIS_DAWN_MODULE_SOURCES})

geant4_module_link_libraries(G4FR
  PUBLIC
    G4geometrymng
    G4modeling
    G4vis_management
    G4globman
  PRIVATE
    G4csg
    G4graphics_reps
    G4hepgeometry)

