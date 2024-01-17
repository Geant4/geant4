# - G4visVTK module build definition
# Define the Geant4 Module
geant4_add_module(G4visVtk
  PUBLIC_HEADERS
    G4Vtk.hh
    G4VtkOffscreen.hh
  PRIVATE_HEADERS
    G4VtkUtility.hh
    G4VtkMessenger.hh
    G4VtkSceneHandler.hh
    G4VtkViewer.hh
    G4VtkOffscreenViewer.hh
    G4VtkInteractorStyle.hh
    G4VtkStore.hh
    G4VtkVisContext.hh
    G4VVtkPipeline.hh
    G4VtkCutterPipeline.hh
    G4VtkClipClosedSurfacePipeline.hh
    G4VtkClipOpenPipeline.hh
    G4VtkImagePipeline.hh
    G4VtkPolydataPipeline.hh
    G4VtkPolydataInstancePipeline.hh
    G4VtkPolydataInstanceTensorPipeline.hh
    G4VtkPolydataInstanceAppendPipeline.hh
    G4VtkPolydataInstanceBakePipeline.hh
    G4VtkPolydataPolylinePipeline.hh
    G4VtkPolydataPolyline2DPipeline.hh
    G4VtkPolydataSpherePipeline.hh
    G4VtkPolydataCubePipeline.hh
    G4VtkPolydataPolygonPipeline.hh
    G4VtkTextPipeline.hh
    G4VtkText2DPipeline.hh
    G4VtkStructuredGridPipeline.hh
    G4VtkUnstructuredGridPipeline.hh
    vtkTensorGlyphColor.h

  SOURCES
    G4Vtk.cc
    G4VtkOffscreen.cc
    G4VtkUtility.cc
    G4VtkMessenger.cc
    G4VtkSceneHandler.cc
    G4VtkViewer.cc
    G4VtkOffscreenViewer.cc
    G4VtkInteractorStyle.cc
    G4VtkStore.cc
    G4VtkCutterPipeline.cc
    G4VtkClipClosedSurfacePipeline.cc
    G4VtkClipOpenPipeline.cc
    G4VtkImagePipeline.cc
    G4VtkPolydataPipeline.cc
    G4VtkPolydataInstancePipeline.cc
    G4VtkPolydataInstanceTensorPipeline.cc
    G4VtkPolydataInstanceAppendPipeline.cc
    G4VtkPolydataInstanceBakePipeline.cc
    G4VtkPolydataPolylinePipeline.cc
    G4VtkPolydataPolyline2DPipeline.cc
    G4VtkPolydataSpherePipeline.cc
    G4VtkPolydataPolygonPipeline.cc
    G4VtkTextPipeline.cc
    G4VtkText2DPipeline.cc
    G4VtkUnstructuredGridPipeline.cc
    vtkTensorGlyphColor.cxx)

geant4_module_compile_definitions(G4visVtk PUBLIC G4VIS_USE_VTK)

geant4_module_link_libraries(G4visVtk
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
    ${VTK_LIBRARIES})

# - VTK-Qt if Qt enabled
if(GEANT4_USE_QT)
  geant4_module_sources(G4visVtk
    PUBLIC_HEADERS
      G4VtkQt.hh
    PRIVATE_HEADERS
      G4VtkQtSceneHandler.hh
      G4VtkQtViewer.hh
    SOURCES
      G4VtkQt.cc
      G4VtkQtSceneHandler.cc
      G4VtkQtViewer.cc)

  geant4_module_compile_definitions(G4visVtk PUBLIC G4VIS_USE_VTK_QT)
  geant4_module_link_libraries(G4visVtk PRIVATE G4UIimplementation)
endif()

