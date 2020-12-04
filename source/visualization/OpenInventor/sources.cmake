#------------------------------------------------------------------------------
# Module : G4OpenInventor
# Package: Geant4.src.G4visualization.G4OpenInventor
#------------------------------------------------------------------------------

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
add_definitions(-DG4VIS_BUILD_OI_DRIVER)

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
  add_definitions(-DG4INTY_BUILD_XT)
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

  # Add the definitions
  # Argh.. Have to remember about INTY and UI because of their use...
  add_definitions(-DG4VIS_BUILD_OIQT_DRIVER)
  add_definitions(-DG4INTY_BUILD_QT)
  add_definitions(-DG4UI_BUILD_QT_SESSION)

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

  # Add the definitions for SoWin
  add_definitions(-DG4INTY_BUILD_WIN32)
  add_definitions(-DG4VIS_BUILD_OIWIN32_DRIVER)

  # SoWin Library
  list(APPEND G4VIS_MODULE_OPENINVENTOR_LINK_LIBRARIES SoWin::SoWin OpenGL::GL)
endif()

#
# Define the Geant4 Module.
#
geant4_define_module(NAME G4OpenInventor
  HEADERS
    ${G4VIS_MODULE_OPENINVENTOR_HEADERS}
  SOURCES
    ${G4VIS_MODULE_OPENINVENTOR_SOURCES}
  GRANULAR_DEPENDENCIES
    G4UIcommon
    G4csg
    G4geometrymng
    G4gl2ps
    G4globman
    G4graphics_reps
    G4hits
    G4intercoms
    G4materials
    G4modeling
    G4specsolids
    G4tracking
    G4vis_management
  GLOBAL_DEPENDENCIES
    G4digits_hits
    G4geometry
    G4gl2ps
    G4global
    G4graphics_reps
    G4intercoms
    G4interfaces
    G4materials
    G4modeling
    G4tracking
    G4vis_management
  LINK_LIBRARIES
    ${G4VIS_MODULE_OPENINVENTOR_LINK_LIBRARIES}
  )

# List any source specific properties here

