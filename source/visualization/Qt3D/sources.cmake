# - G4visQt3D module build definition

# Define the Geant4 Module.
geant4_add_module(G4visQt3D
  PUBLIC_HEADERS
    G4Qt3D.hh
    G4Qt3DQEntity.hh
    G4Qt3DSceneHandler.hh
    G4Qt3DUtils.hh
    G4Qt3DViewer.hh
  SOURCES
    G4Qt3D.cc
    G4Qt3DSceneHandler.cc
    G4Qt3DUtils.cc
    G4Qt3DViewer.cc)

geant4_module_compile_definitions(G4visQt3D PRIVATE G4VIS_BUILD_QT3D_DRIVER)

geant4_module_link_libraries(G4visQt3D
  PUBLIC
    G4modeling
    G4globman
    G4hepgeometry
    G4vis_management
    Qt${QT_VERSION_MAJOR}::Gui
    Qt${QT_VERSION_MAJOR}::Widgets
    Qt${QT_VERSION_MAJOR}::3DCore
    Qt${QT_VERSION_MAJOR}::3DExtras
    Qt${QT_VERSION_MAJOR}::3DRender
  PRIVATE
    G4csg
    G4graphics_reps
    G4geometrymng
    G4navigation
    G4UIbasic
    G4intercoms)




