# - G4visVTK module build definition
# Define the Geant4 Module
geant4_add_module(G4visVtk
  PUBLIC_HEADERS
    G4Vtk.hh
    G4VtkMessenger.hh
    G4VtkSceneHandler.hh
    G4VtkViewer.hh
    vtkTensorGlyphColor.h
  SOURCES
    G4Vtk.cc
    G4VtkMessenger.cc
    G4VtkSceneHandler.cc
    G4VtkViewer.cc
    vtkTensorGlyphColor.cxx)

geant4_module_link_libraries(G4visVtk
  PUBLIC
    G4intercoms
    G4vis_management
    ${VTK_LIBRARIES}
  PRIVATE
    G4csg
    G4geometrymng
    G4globman
    G4graphics_reps
    G4materials
    G4modeling)

# - VTK-Qt if Qt enabled
if(GEANT4_USE_QT)
  geant4_module_sources(G4visVtk
    PUBLIC_HEADERS
      G4VtkQt.hh
      G4VtkQtSceneHandler.hh
      G4VtkQtViewer.hh
    SOURCES
      G4VtkQt.cc
      G4VtkQtSceneHandler.cc
      G4VtkQtViewer.cc)

  geant4_module_link_libraries(G4visVtk PRIVATE G4UIbasic G4UIcommon)
endif()

