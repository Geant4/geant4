# - G4ToolsSG module build definition
# Define the Geant4 Module.
geant4_add_module(G4ToolsSG
  PUBLIC_HEADERS 
    G4ToolsSGOffscreen.hh
  PRIVATE_HEADERS
    G4ToolsSGNode.hh
    G4ToolsSGSceneHandler.hh
    G4ToolsSGViewer.hh
    G4ToolsSGOffscreenViewer.hh
  SOURCES
    G4ToolsSGSceneHandler.cc
    G4ToolsSGOffscreen.cc)

geant4_module_link_libraries(G4ToolsSG
  PUBLIC
    G4vis_management
  PRIVATE
    G4graphics_reps
    G4intercoms
    G4modeling
    G4navigation
    G4tools)

# Freetype support if selected
if(GEANT4_USE_FREETYPE)
  geant4_module_compile_definitions(G4ToolsSG PRIVATE TOOLS_USE_FREETYPE)
  geant4_module_link_libraries(G4ToolsSG PUBLIC Freetype::Freetype)
endif()

# X11 sources if selected
if(GEANT4_USE_TOOLSSG_X11_GLES)
  geant4_module_sources(G4ToolsSG PUBLIC_HEADERS G4ToolsSGX11GLES.hh SOURCES G4ToolsSGX11GLES.cc)
  geant4_module_compile_definitions(G4ToolsSG PUBLIC G4VIS_USE_TOOLSSG_X11_GLES)
endif()

if(GEANT4_USE_TOOLSSG_X11_ZB)
  geant4_module_sources(G4ToolsSG PUBLIC_HEADERS G4ToolsSGX11ZB.hh SOURCES G4ToolsSGX11ZB.cc)
  geant4_module_compile_definitions(G4ToolsSG PUBLIC G4VIS_USE_TOOLSSG_X11_ZB)
endif()

# Xt sources if selected
if(GEANT4_USE_TOOLSSG_XT_GLES)
  geant4_module_sources(G4ToolsSG PUBLIC_HEADERS G4ToolsSGXtGLES.hh SOURCES G4ToolsSGXtGLES.cc)
  geant4_module_compile_definitions(G4ToolsSG 
    PUBLIC G4VIS_USE_TOOLSSG_XT_GLES
    PRIVATE TOOLS_USE_GL_GL_H)
  geant4_module_link_libraries(G4ToolsSG PRIVATE G4UIimplementation)
endif()

if(GEANT4_USE_TOOLSSG_XT_ZB)
  geant4_module_sources(G4ToolsSG PUBLIC_HEADERS G4ToolsSGXtZB.hh SOURCES G4ToolsSGXtZB.cc)
  geant4_module_compile_definitions(G4ToolsSG PUBLIC G4VIS_USE_TOOLSSG_XT_ZB)
  geant4_module_link_libraries(G4ToolsSG PRIVATE G4UIimplementation)
endif()

# X11/Xt links if selected
if(GEANT4_USE_TOOLSSG_X11_GLES OR GEANT4_USE_TOOLSSG_XT_GLES)
  geant4_module_link_libraries(G4ToolsSG PRIVATE X11::SM X11::ICE X11::X11 X11::Xext X11::Xmu)
  if(APPLE)
    geant4_module_link_libraries(G4ToolsSG PRIVATE XQuartzGL::GL)
  else()
    geant4_module_link_libraries(G4ToolsSG PRIVATE OpenGL::GL)
  endif()
endif()

if(GEANT4_USE_TOOLSSG_X11_ZB OR GEANT4_USE_TOOLSSG_XT_ZB)
  geant4_module_link_libraries(G4ToolsSG PRIVATE X11::SM X11::ICE X11::X11 X11::Xext X11::Xmu)
endif()

# Qt sources/links if selected
if(GEANT4_USE_TOOLSSG_QT_GLES)
  geant4_module_sources(G4ToolsSG
    PUBLIC_HEADERS 
      G4ToolsSGQtGLES.hh
    PRIVATE_HEADERS
      G4ToolsSGQtGLESViewer.hh
    SOURCES
      G4ToolsSGQtGLES.cc
      G4ToolsSGQtGLESViewer.cc)

  geant4_module_compile_definitions(G4ToolsSG PUBLIC G4VIS_USE_TOOLSSG_QT_GLES)

  geant4_module_link_libraries(G4ToolsSG 
    PRIVATE 
      Qt${QT_VERSION_MAJOR}::OpenGL
      Qt${QT_VERSION_MAJOR}::Gui
      Qt${QT_VERSION_MAJOR}::Widgets 
      OpenGL::GL
      G4UIimplementation)
  
  if(QT_VERSION_MAJOR GREATER 5)
    geant4_module_link_libraries(G4ToolsSG PRIVATE Qt${QT_VERSION_MAJOR}::OpenGLWidgets)
  endif()

  geant4_set_module_property(G4ToolsSG PROPERTY AUTOMOC ON)

  # Minor hack for MOC-ing. Qt's moc requires visibility of the private headers
  # - Will not affect external consumers and should be minimal impact interanally
  #   as this is a leaf category
  geant4_module_include_directories(G4ToolsSG
    PRIVATE 
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/private>)
endif()

if(GEANT4_USE_TOOLSSG_QT_ZB)
  geant4_module_sources(G4ToolsSG
    PUBLIC_HEADERS
      G4ToolsSGQtZB.hh
    PRIVATE_HEADERS
      G4ToolsSGQtZBViewer.hh
    SOURCES
      G4ToolsSGQtZB.cc
      G4ToolsSGQtZBViewer.cc)

  geant4_module_compile_definitions(G4ToolsSG PUBLIC G4VIS_USE_TOOLSSG_QT_ZB)

  geant4_module_link_libraries(G4ToolsSG
    PRIVATE
      Qt${QT_VERSION_MAJOR}::Gui
      Qt${QT_VERSION_MAJOR}::Widgets
      G4UIimplementation)

  geant4_set_module_property(G4ToolsSG PROPERTY AUTOMOC ON)

  # Minor hack for MOC-ing. Qt's moc requires visibility of the private headers
  # - Will not affect external consumers and should be minimal impact interanally
  #   as this is a leaf category
  geant4_module_include_directories(G4ToolsSG
    PRIVATE
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/private>)
endif()

# Windows sources if selected
if(GEANT4_USE_TOOLSSG_WINDOWS_GLES)
  geant4_module_sources(G4ToolsSG PUBLIC_HEADERS G4ToolsSGWindowsGLES.hh SOURCES G4ToolsSGWindowsGLES.cc)
  geant4_module_compile_definitions(G4ToolsSG PUBLIC G4VIS_USE_TOOLSSG_WINDOWS_GLES)
  geant4_module_link_libraries(G4ToolsSG PRIVATE OpenGL::GL)
endif()

if(GEANT4_USE_TOOLSSG_WINDOWS_ZB)
  geant4_module_sources(G4ToolsSG PUBLIC_HEADERS G4ToolsSGWindowsZB.hh SOURCES G4ToolsSGWindowsZB.cc)
  geant4_module_compile_definitions(G4ToolsSG PUBLIC G4VIS_USE_TOOLSSG_WINDOWS_ZB)
endif()

