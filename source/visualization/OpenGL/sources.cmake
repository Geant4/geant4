# - G4OpenGL module build definition

# Define the Geant4 Module.
geant4_add_module(G4OpenGL
  PUBLIC_HEADERS 
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
    G4gl2ps.hh
  SOURCES
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
    G4gl2ps.cc ) 

geant4_module_link_libraries(G4OpenGL
  PUBLIC
    G4globman
    G4modeling
    G4graphics_reps
    G4hepgeometry
    G4intercoms
    G4vis_management
  PRIVATE
    G4geometrymng
    G4run)

#----------------------------------------------------------------------------
# Add X11 OpenGL Support if requested
if(GEANT4_USE_OPENGL_X11)
  geant4_module_sources(G4OpenGL
    PUBLIC_HEADERS
      G4OpenGLImmediateX.hh
      G4OpenGLImmediateXViewer.hh
      G4OpenGLStoredX.hh
      G4OpenGLStoredXViewer.hh
      G4OpenGLXViewer.hh
    SOURCES
      G4OpenGLImmediateX.cc
      G4OpenGLImmediateXViewer.cc
      G4OpenGLStoredX.cc
      G4OpenGLStoredXViewer.cc
      G4OpenGLXViewer.cc)

  # Add the compile definitions needed for the X11 component (G4OpenGL.hh, G4OpenGLViewer.cc)
  geant4_module_compile_definitions(G4OpenGL PRIVATE G4VIS_BUILD_OPENGLX_DRIVER)
endif()

#----------------------------------------------------------------------------
# Add Xm OpenGL Support if requested
if(GEANT4_USE_XM)
  geant4_module_sources(G4OpenGL
    PUBLIC_HEADERS
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
    SOURCES
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
      G4OpenGLXmWindowHandlingCallbacks.cc )

  # Special case of building Xm without X11
  if(NOT GEANT4_USE_OPENGL_X11)
    geant4_module_sources(G4OpenGL PUBLIC_HEADERS G4OpenGLXViewer.hh SOURCES G4OpenGLXViewer.cc)
    geant4_module_compile_definitions(G4OpenGL PRIVATE G4VIS_BUILD_OPENGLX_DRIVER)
  endif()

  # Add the compile definitions needed for the Xm component (G4OpenGL.hh, G4OpenGLViewer.cc)
  geant4_module_compile_definitions(G4OpenGL PRIVATE G4VIS_BUILD_OPENGLXM_DRIVER)

  # Add in Xm and needed modules
  geant4_module_link_libraries(G4OpenGL PUBLIC Motif::Xm PRIVATE G4UIcommon)
endif()

# Common X11/Xm link libraries
if(GEANT4_USE_OPENGL_X11 OR GEANT4_USE_XM)
  geant4_module_link_libraries(G4OpenGL PUBLIC X11::SM X11::ICE X11::X11 X11::Xext X11::Xmu)
  if(APPLE)
    geant4_module_link_libraries(G4OpenGL PUBLIC XQuartzGL::GL)
  else()
    geant4_module_link_libraries(G4OpenGL PUBLIC OpenGL::GL)
  endif() 
endif()  

#----------------------------------------------------------------------------
# Add Qt OpenGL support if requested
if(GEANT4_USE_QT)
  geant4_module_sources(G4OpenGL
    PUBLIC_HEADERS
      G4OpenGLImmediateQt.hh
      G4OpenGLImmediateQtViewer.hh
      G4OpenGLQt.hh
      G4OpenGLQtExportDialog.hh
      G4OpenGLQtMovieDialog.hh
      G4OpenGLVboDrawer.hh
      G4OpenGLQtViewer.hh
      G4OpenGLStoredQt.hh
      G4OpenGLStoredQtSceneHandler.hh
      G4OpenGLStoredQtViewer.hh
    SOURCES
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

  # Add the definitions (G4OpenGL.hh, G4OpenGLViewer.cc)
  geant4_module_compile_definitions(G4OpenGL PRIVATE G4VIS_BUILD_OPENGLQT_DRIVER)

  # Add in Qt libraries and geant4 modules
  geant4_module_link_libraries(G4OpenGL
    PUBLIC Qt${QT_VERSION_MAJOR}::OpenGL Qt${QT_VERSION_MAJOR}::Gui Qt${QT_VERSION_MAJOR}::PrintSupport Qt${QT_VERSION_MAJOR}::Widgets OpenGL::GL
    PRIVATE G4UIbasic G4UIcommon)
endif()

#----------------------------------------------------------------------------
# Add WIN32 support, if requested
#
if(GEANT4_USE_OPENGL_WIN32)
  geant4_module_sources(G4OpenGL
    PUBLIC_HEADERS
      G4OpenGLImmediateWin32.hh
      G4OpenGLImmediateWin32Viewer.hh
      G4OpenGLStoredWin32.hh
      G4OpenGLStoredWin32Viewer.hh
      G4OpenGLWin32Viewer.hh
    SOURCES
      G4OpenGLImmediateWin32.cc
      G4OpenGLImmediateWin32Viewer.cc
      G4OpenGLStoredWin32.cc
      G4OpenGLStoredWin32Viewer.cc
      G4OpenGLWin32Viewer.cc )

  # Add the compile definitions (G4OpenGL.hh, G4OpenGLViewer.cc)
  geant4_module_compile_definitions(G4OpenGL PRIVATE G4VIS_BUILD_OPENGLWIN32_DRIVER)
  geant4_module_link_libraries(G4OpenGL PUBLIC OpenGL::GL)
endif()

