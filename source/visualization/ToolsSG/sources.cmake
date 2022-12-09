# - G4ToolsSG module build definition
# Define the Geant4 Module.
geant4_add_module(G4ToolsSG
  PUBLIC_HEADERS 
    G4ToolsSGNode.hh
    G4ToolsSGSceneHandler.hh
    G4ToolsSGViewer.hh
    G4ToolsSGOffscreen.hh
    G4ToolsSGOffscreenViewer.hh
  SOURCES
    G4ToolsSGSceneHandler.cc
    G4ToolsSGOffscreen.cc)

geant4_module_link_libraries(G4ToolsSG
  PUBLIC
    G4intercoms
    G4modeling
    G4vis_management
    G4tools
  PRIVATE
    G4graphics_reps
    G4navigation)

# Freetype support if selected
if(GEANT4_USE_FREETYPE)
  geant4_module_compile_definitions(G4ToolsSG PRIVATE TOOLS_USE_FREETYPE)
  geant4_module_link_libraries(G4ToolsSG PUBLIC Freetype::Freetype)
endif()

# X11 sources if selected
if(GEANT4_USE_TOOLSSG_X11_GLES)
  geant4_module_sources(G4ToolsSG PUBLIC_HEADERS G4ToolsSGX11GLES.hh SOURCES G4ToolsSGX11GLES.cc)
endif()

# Xt sources if selected
if(GEANT4_USE_TOOLSSG_XT_GLES)
  geant4_module_sources(G4ToolsSG PUBLIC_HEADERS G4ToolsSGXtGLES.hh SOURCES G4ToolsSGXtGLES.cc)
  geant4_module_compile_definitions(G4ToolsSG PRIVATE TOOLS_USE_GL_GL_H)

  # A minor hack around a likely issue in geant4_module_link_libraries
  # When Qt is activated, G4UIcommon is a public dependency. CMake will deduplicate
  # this (in favour of PUBLIC at final library link time, but geant4_module_link_libraries
  # does not do this internally yet (To be checked). We then get validation
  # warnings on duplicated deps. NB: Also demonstrates awkward vis system of more than
  # one driver per library...
  if(GEANT4_USE_TOOLSSG_QT_GLES)
    geant4_module_link_libraries(G4ToolsSG PUBLIC G4UIcommon)
  else()
    geant4_module_link_libraries(G4ToolsSG PRIVATE G4UIcommon)
  endif()
endif()

# X11/Xt links if selected
if(GEANT4_USE_TOOLSSG_X11_GLES OR GEANT4_USE_TOOLSSG_XT_GLES)
  geant4_module_link_libraries(G4ToolsSG PUBLIC X11::SM X11::ICE X11::X11 X11::Xext X11::Xmu)
  if(APPLE)
    geant4_module_link_libraries(G4ToolsSG PUBLIC XQuartzGL::GL)
  else()
    geant4_module_link_libraries(G4ToolsSG PUBLIC OpenGL::GL)
  endif()
endif()

# Qt sources/links if selected
if(GEANT4_USE_TOOLSSG_QT_GLES)
  geant4_module_sources(G4ToolsSG
    PUBLIC_HEADERS 
      G4ToolsSGQtViewer.hh
      G4ToolsSGQtGLES.hh
    SOURCES
      G4ToolsSGQtGLES.cc
      G4ToolsSGQtViewer.cc)

  geant4_module_link_libraries(G4ToolsSG 
    PUBLIC 
      Qt${QT_VERSION_MAJOR}::OpenGL
      Qt${QT_VERSION_MAJOR}::Gui
      Qt${QT_VERSION_MAJOR}::PrintSupport
      Qt${QT_VERSION_MAJOR}::Widgets 
      OpenGL::GL
      G4UIbasic
      G4UIcommon)
endif()

# Windows sources if selected
if(GEANT4_USE_TOOLSSG_WINDOWS_GLES)
  geant4_module_sources(G4ToolsSG PUBLIC_HEADERS G4ToolsSGWindowsGLES.hh SOURCES G4ToolsSGWindowsGLES.cc)
  geant4_module_link_libraries(G4ToolsSG PUBLIC OpenGL::GL)
endif()




