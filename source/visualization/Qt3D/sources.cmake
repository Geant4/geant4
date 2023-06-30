# - G4visQt3D module build definition

# Define the Geant4 Module.
geant4_add_module(G4visQt3D
  PUBLIC_HEADERS
    G4Qt3D.hh
  PRIVATE_HEADERS
    G4Qt3DQEntity.hh
    G4Qt3DSceneHandler.hh
    G4Qt3DUtils.hh
    G4Qt3DViewer.hh
  SOURCES
    G4Qt3D.cc
    G4Qt3DSceneHandler.cc
    G4Qt3DUtils.cc
    G4Qt3DViewer.cc)

geant4_module_compile_definitions(G4visQt3D PUBLIC G4VIS_USE_QT3D)

geant4_module_link_libraries(G4visQt3D
  PUBLIC
    G4vis_management
  PRIVATE
    G4csg
    G4graphics_reps
    G4globman
    G4geometrymng
    G4hepgeometry
    G4modeling
    G4navigation
    G4UIimplementation
    G4intercoms
    Qt${QT_VERSION_MAJOR}::Gui
    Qt${QT_VERSION_MAJOR}::Widgets
    Qt${QT_VERSION_MAJOR}::3DCore
    Qt${QT_VERSION_MAJOR}::3DExtras
    Qt${QT_VERSION_MAJOR}::3DRender)
