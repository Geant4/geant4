#------------------------------------------------------------------------------
# Module : G4FR
# Package: Geant4.src.G4visualization.G4FR
#------------------------------------------------------------------------------

#
# Module has optional sources
#

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

#
# Define the Geant4 Module.
#
geant4_define_module(NAME G4FR
  HEADERS
    ${G4VIS_DAWN_MODULE_HEADERS}
  SOURCES
    ${G4VIS_DAWN_MODULE_SOURCES}
  GRANULAR_DEPENDENCIES
    G4csg
    G4geometrymng
    G4globman
    G4graphics_reps
    G4hits
    G4intercoms
    G4modeling
    G4specsolids
    G4vis_management
  GLOBAL_DEPENDENCIES
    G4digits_hits
    G4geometry
    G4global
    G4graphics_reps
    G4intercoms
    G4modeling
    G4vis_management
)

# List any source specific properties here
