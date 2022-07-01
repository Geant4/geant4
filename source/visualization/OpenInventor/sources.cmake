# - G4OpenInventor module build definition
# Define the Geant4 Module.
geant4_add_module(G4OpenInventor
  PUBLIC_HEADERS
    G4OpenInventor.hh
    G4OpenInventorSceneHandler.hh
    G4OpenInventorTransform3D.hh
    G4OpenInventorViewer.hh
    G4VisFeaturesOfOpenInventor.hh
    Geant4_SoPolyhedron.h
    SoG4LineSet.h
    SoG4MarkerSet.h
    SoG4Polyhedron.h
  SOURCES
    G4OpenInventor.cc
    G4OpenInventorSceneHandler.cc
    G4OpenInventorTransform3D.cc
    G4OpenInventorViewer.cc
    G4VisFeaturesOfOpenInventor.cc
    SbPainter.cc
    SbPainterPS.cc
    SoAlternateRepAction.cc
    SoBox.cc
    SoCons.cc
    SoCounterAction.cc
    SoDetectorTreeKit.cc
    SoGL2PSAction.cc
    SoImageWriter.cc
    SoMarkerSet.cc
    SoPolyhedron.cc
    SoStyleCache.cc
    SoTrap.cc
    SoTrd.cc
    SoTubs.cc )

geant4_module_link_libraries(G4OpenInventor
  PUBLIC
    G4UIcommon
    G4globman
    G4graphics_reps
    G4hepgeometry
    G4intercoms
    G4modeling
    G4vis_management
    Coin::Coin
  PRIVATE
    G4UIbasic
    G4csg
    G4geometrymng
    G4tools
    G4materials
    G4tracking)

# UNIX Only (Xt) sources
if(GEANT4_USE_INVENTOR_XT)
  geant4_module_sources(G4OpenInventor 
    PUBLIC_HEADERS
      G4OpenInventorX.hh
      G4OpenInventorXt.hh
      G4OpenInventorXtExaminerViewerMessenger.hh
      G4OpenInventorXtExaminerViewer.hh
      G4OpenInventorXtExtended.hh
      G4OpenInventorXtExtendedViewer.hh
      G4OpenInventorXtViewer.hh
      wheelmouse.h
      SoXtInternal.h
      console.h
      favorites.h
      saveViewPt.h
      pickext.h
      pickref.h
      wireframe.h
    SOURCES
      G4OpenInventorXt.cc
      G4OpenInventorXtExaminerViewer.cc
      G4OpenInventorXtExaminerViewerMessenger.cc
      G4OpenInventorXtExtended.cc
      G4OpenInventorXtExtendedViewer.cc
      G4OpenInventorXtViewer.cc
      wheelmouse.cc )
  
  geant4_module_compile_definitions(G4OpenInventor PRIVATE G4VIS_BUILD_OIX_DRIVER)
  geant4_module_link_libraries(G4OpenInventor PUBLIC SoXt::SoXt Motif::Xm)
  if(APPLE)
    geant4_module_link_libraries(G4OpenInventor PUBLIC XQuartzGL::GL)
  else()
    geant4_module_link_libraries(G4OpenInventor PUBLIC OpenGL::GL)
  endif()
endif()

# Qt-only sources
if(GEANT4_USE_INVENTOR_QT)
  geant4_module_sources(G4OpenInventor
    PUBLIC_HEADERS
      G4OpenInventorQt.hh
      G4OpenInventorQtExaminerViewer.hh
      G4OpenInventorQtViewer.hh
      G4SoQt.hh
      ui_OIQtListsDialog.h
    SOURCES
      G4OpenInventorQt.cc
      G4OpenInventorQtExaminerViewer.cc
      G4OpenInventorQtViewer.cc
      G4SoQt.cc )

  geant4_module_link_libraries(G4OpenInventor PUBLIC SoQt::SoQt OpenGL::GL Qt${QT_VERSION_MAJOR}::OpenGL Qt${QT_VERSION_MAJOR}::Gui Qt${QT_VERSION_MAJOR}::PrintSupport Qt${QT_VERSION_MAJOR}::Widgets)
endif()

# - WIN32 Only (Win32) sources
if(GEANT4_USE_INVENTOR_WIN)
  geant4_module_sources(G4OpenInventor
    PUBLIC_HEADERS
      G4OpenInventorWin.hh
      G4OpenInventorWin32.hh
      G4OpenInventorWinViewer.hh
    SOURCES
      G4OpenInventorWin.cc
      G4OpenInventorWinViewer.cc)
  
  geant4_module_link_libraries(G4OpenInventor PUBLIC SoWin::SoWin OpenGL::GL)
endif()
