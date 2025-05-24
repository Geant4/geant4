# - G4visVTK module build definition
# Define the Geant4 Module
geant4_add_module(G4visXr
  PUBLIC_HEADERS
    G4Xr.hh
    G4XrQt.hh
  PRIVATE_HEADERS
    G4XrViewer.hh
    G4XrQtViewer.hh
    G4XrSceneHandler.hh
  SOURCES
    G4Xr.cc
    G4XrQt.cc
    G4XrViewer.cc
    G4XrQtViewer.cc
    G4XrSceneHandler.cc)

geant4_module_compile_definitions(G4visXr PUBLIC G4VIS_USE_XR)

geant4_module_link_libraries(G4visXr
  PUBLIC
    G4vis_management
  PRIVATE
    G4csg
    G4geometrymng
    G4globman
    G4graphics_reps
    G4hepgeometry
    G4intercoms
    G4materials
    G4modeling
    G4specsolids
    Qt${QT_VERSION_MAJOR}::Gui
    Qt${QT_VERSION_MAJOR}::Widgets
    Qt${QT_VERSION_MAJOR}::3DCore
    Qt${QT_VERSION_MAJOR}::3DExtras
    Qt${QT_VERSION_MAJOR}::3DRender
    G4tinygltf
    G4cpp-httplib)

geant4_module_link_libraries(G4visXr PRIVATE G4UIimplementation)

