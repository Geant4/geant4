# - G4OpenInventor module build definition

#----------------------------------------------------------------------------
# Geant4 OpenInventor Core sources and headers (all platforms)
#
set(G4VIS_MODULE_OPENINVENTOR_HEADERS
  G4OpenInventor.hh
  G4OpenInventorSceneHandler.hh
  G4OpenInventorTransform3D.hh
  G4OpenInventorViewer.hh
  G4VisFeaturesOfOpenInventor.hh
  Geant4_SoPolyhedron.h
  SoG4LineSet.h
  SoG4MarkerSet.h
  SoG4Polyhedron.h
  )

set(G4VIS_MODULE_OPENINVENTOR_SOURCES
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
  SoTubs.cc
  )

set(G4VIS_MODULE_OPENINVENTOR_LINK_LIBRARIES Coin::Coin)

#----------------------------------------------------------------------------
# UNIX Only (Xt) sources
#
if(GEANT4_USE_INVENTOR_XT)
  list(APPEND G4VIS_MODULE_OPENINVENTOR_HEADERS
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
    wireframe.h)

  list(APPEND G4VIS_MODULE_OPENINVENTOR_SOURCES
    G4OpenInventorXt.cc
    G4OpenInventorXtExaminerViewer.cc
    G4OpenInventorXtExaminerViewerMessenger.cc
    G4OpenInventorXtExtended.cc
    G4OpenInventorXtExtendedViewer.cc
    G4OpenInventorXtViewer.cc
    wheelmouse.cc)

  # Add the definitions for SoXt
  add_definitions(-DG4VIS_BUILD_OIX_DRIVER)

  # SoXt Library and others
  list(APPEND G4VIS_MODULE_OPENINVENTOR_LINK_LIBRARIES SoXt::SoXt Motif::Xm)
  if(APPLE)
    list(APPEND G4VIS_MODULE_OPENINVENTOR_LINK_LIBRARIES XQuartzGL::GL)
  else()
    list(APPEND G4VIS_MODULE_OPENINVENTOR_LINK_LIBRARIES OpenGL::GL)
  endif()
endif()

#----------------------------------------------------------------------------
# Open Inventor Qt
#
if(GEANT4_USE_INVENTOR_QT)
  list(APPEND G4VIS_MODULE_OPENINVENTOR_HEADERS
    G4OpenInventorQt.hh
    G4OpenInventorQtExaminerViewer.hh
    G4OpenInventorQtViewer.hh
    G4SoQt.hh
    ui_OIQtListsDialog.h)

  list(APPEND G4VIS_MODULE_OPENINVENTOR_SOURCES
    G4OpenInventorQt.cc
    G4OpenInventorQtExaminerViewer.cc
    G4OpenInventorQtViewer.cc
    G4SoQt.cc)

  # Add libraries
  list(APPEND G4VIS_MODULE_OPENINVENTOR_LINK_LIBRARIES
    SoQt::SoQt
    OpenGL::GL
    Qt5::OpenGL Qt5::Gui Qt5::PrintSupport Qt5::Widgets)
endif()


#----------------------------------------------------------------------------
# WIN32 Only (Win32) sources
#
if(GEANT4_USE_INVENTOR_WIN)
  set(G4VIS_MODULE_OPENINVENTOR_HEADERS
    ${G4VIS_MODULE_OPENINVENTOR_HEADERS}
    G4OpenInventorWin.hh
    G4OpenInventorWin32.hh
    G4OpenInventorWinViewer.hh
    )

  set(G4VIS_MODULE_OPENINVENTOR_SOURCES
    ${G4VIS_MODULE_OPENINVENTOR_SOURCES}
    G4OpenInventorWin.cc
    G4OpenInventorWinViewer.cc
    )

  # SoWin Library
  list(APPEND G4VIS_MODULE_OPENINVENTOR_LINK_LIBRARIES SoWin::SoWin OpenGL::GL)
endif()


# Define the Geant4 Module.
geant4_add_module(G4OpenInventor
  PUBLIC_HEADERS
    ${G4VIS_MODULE_OPENINVENTOR_HEADERS}
  SOURCES
    ${G4VIS_MODULE_OPENINVENTOR_SOURCES})

geant4_module_link_libraries(G4OpenInventor
  PUBLIC
    G4UIcommon
    G4globman
    G4graphics_reps
    G4hepgeometry
    G4intercoms
    G4modeling
    G4vis_management
    ${G4VIS_MODULE_OPENINVENTOR_LINK_LIBRARIES}
  PRIVATE
    G4UIbasic
    G4csg
    G4geometrymng
    G4gl2ps
    G4materials
    G4tracking)