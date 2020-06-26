#------------------------------------------------------------------------------
# Module : G4VRML
# Package: Geant4.src.G4visualization.G4VRML
#------------------------------------------------------------------------------

#
# Module has optional sources
#

# List those always built
set(G4VIS_VRML_MODULE_HEADERS
  G4VRML1File.hh
  G4VRML1FileSceneHandler.hh
  G4VRML1FileViewer.hh
  G4VRML2File.hh
  G4VRML2FileSceneHandler.hh
  G4VRML2FileViewer.hh)

set(G4VIS_VRML_MODULE_SOURCES
  G4VRML1File.cc
  G4VRML1FileSceneHandler.cc
  G4VRML1FileViewer.cc
  G4VRML1SceneHandlerFunc.icc
  G4VRML2File.cc
  G4VRML2FileSceneHandler.cc
  G4VRML2FileViewer.cc
  G4VRML2SceneHandlerFunc.icc)

# VRML Network drivers only if user selected
if(GEANT4_USE_NETWORKVRML)
  list(APPEND G4VIS_VRML_MODULE_HEADERS
    FRClient.h
    G4FRClient.hh
    G4VRML1.hh
    G4VRML1SceneHandler.hh
    G4VRML1Viewer.hh
    G4VRML2.hh
    G4VRML2SceneHandler.hh
    G4VRML2Viewer.hh
    G4VRMLNetConfig.hh)

  list(APPEND G4VIS_VRML_MODULE_SOURCES
    FRClient.cc
    G4FRClient.cc
    G4VRML1.cc
    G4VRML1SceneHandler.cc
    G4VRML1Viewer.cc
    G4VRML2.cc
    G4VRML2SceneHandler.cc
    G4VRML2Viewer.cc)

  #
  # Add extra needed defs here
  #
  add_definitions(-DG4VIS_BUILD_VRML_DRIVER)
endif()

#
# Define the Geant4 Module.
#
geant4_define_module(NAME G4VRML
  HEADERS
    ${G4VIS_VRML_MODULE_HEADERS}
  SOURCES
    ${G4VIS_VRML_MODULE_SOURCES}
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
