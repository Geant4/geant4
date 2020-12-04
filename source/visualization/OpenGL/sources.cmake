#------------------------------------------------------------------------------
# Module : G4OpenGL
# Package: Geant4.src.G4visualization.G4OpenGL
#------------------------------------------------------------------------------

#
# Module has optional sources
#
# Define the core sources, includes and libraries which all Geant4
# OpenGL implementations use
#
set(G4VIS_MODULE_OPENGL_HEADERS
  G4OpenGL.hh
  G4OpenGLImmediateViewer.hh
  G4OpenGLImmediateSceneHandler.hh
  G4OpenGLViewer.hh
  G4OpenGLStoredViewer.hh
  G4OpenGLStoredSceneHandler.hh
  G4OpenGLFontBaseStore.hh
  G4OpenGLSceneHandler.hh
  G4OpenGLSceneHandler.icc
  G4OpenGLViewerMessenger.hh
  G4OpenGLTransform3D.hh
  G4VisFeaturesOfOpenGL.hh
)

set(G4VIS_MODULE_OPENGL_SOURCES
  G4OpenGLImmediateViewer.cc
  G4OpenGLImmediateSceneHandler.cc
  G4OpenGLViewer.cc
  G4OpenGLStoredViewer.cc
  G4OpenGLStoredSceneHandler.cc
  G4OpenGLFontBaseStore.cc
  G4OpenGLSceneHandler.cc
  G4OpenGLViewerMessenger.cc
  G4OpenGLTransform3D.cc
  G4VisFeaturesOfOpenGL.cc
)

add_definitions(-DG4VIS_BUILD_OPENGL_DRIVER)

#----------------------------------------------------------------------------
# Add X11 OpenGL Support if requested
#
if(GEANT4_USE_OPENGL_X11)
  list(APPEND G4VIS_MODULE_OPENGL_HEADERS
    G4OpenGLImmediateX.hh
    G4OpenGLImmediateXViewer.hh
    G4OpenGLStoredX.hh
    G4OpenGLStoredXViewer.hh
    G4OpenGLXViewer.hh)

  list(APPEND G4VIS_MODULE_OPENGL_SOURCES
    G4OpenGLImmediateX.cc
    G4OpenGLImmediateXViewer.cc
    G4OpenGLStoredX.cc
    G4OpenGLStoredXViewer.cc
    G4OpenGLXViewer.cc)

  # Add the compile definitions needed for the X11 component
  add_definitions(-DG4VIS_BUILD_OPENGLX_DRIVER)

  # Add in X11 libraries, plus Xmu library
  list(APPEND G4VIS_MODULE_OPENGL_LINK_LIBRARIES X11::SM X11::ICE X11::X11 X11::Xext X11::Xmu)
  if(APPLE)
    list(APPEND G4VIS_MODULE_OPENGL_LINK_LIBRARIES XQuartzGL::GL)
  else()
    list(APPEND G4VIS_MODULE_OPENGL_LINK_LIBRARIES OpenGL::GL)
  endif()
endif()


#----------------------------------------------------------------------------
# Add Xm OpenGL Support if requested
#
if(GEANT4_USE_XM)
  list(APPEND G4VIS_MODULE_OPENGL_HEADERS
    G4OpenGLImmediateXm.hh
    G4OpenGLImmediateXmViewer.hh
    G4OpenGLStoredXm.hh
    G4OpenGLStoredXmViewer.hh
    G4OpenGLXm.hh
    G4OpenGLXmBox.hh
    G4OpenGLXmFourArrowButtons.hh
    G4OpenGLXmFramedBox.hh
    G4OpenGLXmPushButton.hh
    G4OpenGLXmRadioButton.hh
    G4OpenGLXmResources.hh
    G4OpenGLXmSeparator.hh
    G4OpenGLXmSliderBar.hh
    G4OpenGLXmTextField.hh
    G4OpenGLXmTopLevelShell.hh
    G4OpenGLXmViewer.hh
    G4OpenGLXmViewerMessenger.hh
    G4OpenGLXmVWidgetComponent.hh
    G4OpenGLXmVWidgetContainer.hh
    G4OpenGLXmVWidgetObject.hh
    G4OpenGLXmVWidgetShell.hh
  )

  list(APPEND G4VIS_MODULE_OPENGL_SOURCES
    G4OpenGLImmediateXm.cc
    G4OpenGLImmediateXmViewer.cc
    G4OpenGLStoredXm.cc
    G4OpenGLStoredXmViewer.cc
    G4OpenGLXm.cc
    G4OpenGLXmBox.cc
    G4OpenGLXmConvenienceRoutines.cc
    G4OpenGLXmFourArrowButtons.cc
    G4OpenGLXmFramedBox.cc
    G4OpenGLXmMainMenubarCallbacks.cc
    G4OpenGLXmPanningCallbacks.cc
    G4OpenGLXmPushButton.cc
    G4OpenGLXmRadioButton.cc
    G4OpenGLXmRotationCallbacks.cc
    G4OpenGLXmSeparator.cc
    G4OpenGLXmSliderBar.cc
    G4OpenGLXmStyleCallbacks.cc
    G4OpenGLXmTextField.cc
    G4OpenGLXmTopLevelShell.cc
    G4OpenGLXmViewer.cc
    G4OpenGLXmViewerMessenger.cc
    G4OpenGLXmVWidgetComponent.cc
    G4OpenGLXmVWidgetContainer.cc
    G4OpenGLXmVWidgetObject.cc
    G4OpenGLXmVWidgetShell.cc
    G4OpenGLXmWindowHandlingCallbacks.cc
  )

  # Special case of building Xm without X11
  if(NOT GEANT4_USE_OPENGL_X11)
    list(APPEND G4VIS_MODULE_OPENGL_HEADERS G4OpenGLXViewer.hh)
    list(APPEND G4VIS_MODULE_OPENGL_SOURCES G4OpenGLXViewer.cc)
    add_definitions(-DG4VIS_BUILD_OPENGLX_DRIVER)
  endif()

  # Add the compile definitions needed for the Xm component, remembering
  # to add those for the UI part as well!!!
  add_definitions(-DG4VIS_BUILD_OPENGLXM_DRIVER -DG4INTY_BUILD_XT -DG4UI_BUILD_XM_SESSION)

  # Add in Xm, X11 libraries, plus Xmu library
  set(G4VIS_MODULE_OPENGL_LINK_LIBRARIES Motif::Xm X11::SM X11::ICE X11::X11 X11::Xext X11::Xmu ${G4VIS_MODULE_OPENGL_LINK_LIBRARIES})
  if(APPLE)
    list(APPEND G4VIS_MODULE_OPENGL_LINK_LIBRARIES XQuartzGL::GL)
  else()
    set(G4VIS_MODULE_OPENGL_LINK_LIBRARIES OpenGL::GL)
  endif()
endif()

#----------------------------------------------------------------------------
# Add Qt OpenGL support if requested
#
if(GEANT4_USE_QT)
  list(APPEND G4VIS_MODULE_OPENGL_HEADERS
    G4OpenGLImmediateQt.hh
    G4OpenGLImmediateQtViewer.hh
    G4OpenGLQt.hh
    G4OpenGLQtExportDialog.hh
    G4OpenGLQtMovieDialog.hh
    G4OpenGLVboDrawer.hh
    G4OpenGLQtViewer.hh
    G4OpenGLStoredQt.hh
    G4OpenGLStoredQtSceneHandler.hh
    G4OpenGLStoredQtViewer.hh)

  list(APPEND G4VIS_MODULE_OPENGL_SOURCES
    G4OpenGLImmediateQt.cc
    G4OpenGLImmediateQtViewer.cc
    G4OpenGLQt.cc
    G4OpenGLQtExportDialog.cc
    G4OpenGLQtMovieDialog.cc
    G4OpenGLVboDrawer.cc
    G4OpenGLQtViewer.cc
    G4OpenGLStoredQt.cc
    G4OpenGLStoredQtSceneHandler.cc
    G4OpenGLStoredQtViewer.cc)

  # Add the definitions
  # Argh.. Have to remember about INTY and UI because of their use...
  add_definitions(-DG4VIS_BUILD_OPENGLQT_DRIVER -DG4INTY_BUILD_QT -DG4UI_BUILD_QT_SESSION)

  # Add in Qt libraries
  list(APPEND G4VIS_MODULE_OPENGL_LINK_LIBRARIES Qt5::OpenGL Qt5::Gui Qt5::PrintSupport Qt5::Widgets OpenGL::GL)
endif()


#----------------------------------------------------------------------------
# Add WIN32 support, if requested
#
if(GEANT4_USE_OPENGL_WIN32)
  list(APPEND G4VIS_MODULE_OPENGL_HEADERS
    G4OpenGLImmediateWin32.hh
    G4OpenGLImmediateWin32Viewer.hh
    G4OpenGLStoredWin32.hh
    G4OpenGLStoredWin32Viewer.hh
    G4OpenGLWin32Viewer.hh
  )

  list(APPEND G4VIS_MODULE_OPENGL_SOURCES
    G4OpenGLImmediateWin32.cc
    G4OpenGLImmediateWin32Viewer.cc
    G4OpenGLStoredWin32.cc
    G4OpenGLStoredWin32Viewer.cc
    G4OpenGLWin32Viewer.cc
  )

  # Add the compile definitions
  add_definitions(-DG4VIS_BUILD_OPENGLWIN32_DRIVER)
  list(APPEND G4VIS_MODULE_OPENGL_LINK_LIBRARIES OpenGL::GL)
endif()

# May have duplicates in link list, so remove
list(REMOVE_DUPLICATES G4VIS_MODULE_OPENGL_LINK_LIBRARIES)

#----------------------------------------------------------------------------
# Define the Geant4 Module.
#
geant4_define_module(NAME G4OpenGL
  HEADERS
    ${G4VIS_MODULE_OPENGL_HEADERS}
  SOURCES
    ${G4VIS_MODULE_OPENGL_SOURCES}
  GRANULAR_DEPENDENCIES
    G4UIbasic
    G4UIcommon
    G4csg
    G4event
    G4run
    G4tasking
    G4particles
    G4processes
    G4track
    G4materials
    G4geometrymng
    G4gl2ps
    G4globman
    G4graphics_reps
    G4digits
    G4hits
    G4intercoms
    G4modeling
    G4specsolids
    G4tracking
    G4vis_management
  GLOBAL_DEPENDENCIES
    G4digits_hits
    G4event
    G4run
    G4tasking
    G4particles
    G4processes
    G4track
    G4materials
    G4geometry
    G4gl2ps
    G4global
    G4graphics_reps
    G4intercoms
    G4interfaces
    G4modeling
    G4tracking
    G4vis_management
  LINK_LIBRARIES
    ${G4VIS_MODULE_OPENGL_LINK_LIBRARIES}
)

# List any source specific properties here

