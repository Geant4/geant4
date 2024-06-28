# - G4OpenInventor module build definition
# Define the Geant4 Module.

#include "G4OpenInventorX.hh"          // no_geant4_module_check
#include "G4OpenInventorXtExtended.hh" // no_geant4_module_check
#include "G4OpenInventorQt.hh" // no_geant4_module_check
#include "G4OpenInventorWin32.hh" // no_geant4_module_check

geant4_add_module(G4OpenInventor
  PUBLIC_HEADERS
    G4OpenInventor.hh
  PRIVATE_HEADERS
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

# Needed by G4VisExecutive...
geant4_module_compile_definitions(G4OpenInventor PUBLIC G4VIS_USE_OI)

geant4_module_link_libraries(G4OpenInventor
  PUBLIC
    G4vis_management
  PRIVATE
    G4UIcore
    G4UIimplementation
    G4csg
    G4geometrymng
    G4globman
    G4graphics_reps
    G4hepgeometry
    G4intercoms
    G4materials
    G4modeling
    G4tools
    G4tracking
    Coin::Coin)

# UNIX Only (Xt) sources
if(GEANT4_USE_INVENTOR_XT)
  geant4_module_sources(G4OpenInventor 
    PUBLIC_HEADERS
      G4OpenInventorX.hh
      G4OpenInventorXt.hh
      G4OpenInventorXtExtended.hh
    PRIVATE_HEADERS
      G4OpenInventorXtExaminerViewerMessenger.hh
      G4OpenInventorXtExaminerViewer.hh
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
  
  geant4_module_compile_definitions(G4OpenInventor
    PUBLIC G4VIS_USE_OIX
    PRIVATE G4VIS_BUILD_OIX_DRIVER)

  geant4_module_link_libraries(G4OpenInventor PRIVATE SoXt::SoXt Motif::Xm)
  if(APPLE)
    geant4_module_link_libraries(G4OpenInventor PRIVATE XQuartzGL::GL)
  else()
    geant4_module_link_libraries(G4OpenInventor PRIVATE OpenGL::GL)
  endif()
endif()

# Qt-only sources
if(GEANT4_USE_INVENTOR_QT)
  geant4_module_sources(G4OpenInventor
    PUBLIC_HEADERS
      G4OpenInventorQt.hh
    PRIVATE_HEADERS
      G4OpenInventorQtExaminerViewer.hh
      G4OpenInventorQtViewer.hh
      G4SoQt.hh
      ui_OIQtListsDialog.h
    SOURCES
      G4OpenInventorQt.cc
      G4OpenInventorQtExaminerViewer.cc
      G4OpenInventorQtViewer.cc
      G4SoQt.cc )

  geant4_module_compile_definitions(G4OpenInventor PUBLIC G4VIS_USE_OIQT)

  geant4_module_link_libraries(G4OpenInventor PRIVATE SoQt::SoQt OpenGL::GL Qt${QT_VERSION_MAJOR}::OpenGL Qt${QT_VERSION_MAJOR}::Gui Qt${QT_VERSION_MAJOR}::Widgets)

  geant4_set_module_property(G4OpenInventor PROPERTY AUTOMOC ON)
  # Minor hack for MOC-ing. Qt's moc requires visibility of the private headers
  # - Will not affect external consumers and should be minimal impact interanally
  #   as this is a leaf category
  geant4_module_include_directories(G4OpenInventor
    PRIVATE
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/private>)
endif()

# - WIN32 Only (Win32) sources
if(GEANT4_USE_INVENTOR_WIN)
  geant4_module_sources(G4OpenInventor
    PUBLIC_HEADERS
      G4OpenInventorWin.hh
      G4OpenInventorWin32.hh
    PRIVATE_HEADERS
      G4OpenInventorWinViewer.hh
    SOURCES
      G4OpenInventorWin.cc
      G4OpenInventorWinViewer.cc)
  
  geant4_module_compile_definitions(G4OpenInventor PUBLIC G4VIS_USE_OIWIN32)
  
  geant4_module_link_libraries(G4OpenInventor PRIVATE SoWin::SoWin OpenGL::GL)
endif()
