# - G4ToolsSG module build definition
# Define the Geant4 Module.
geant4_add_module(G4ToolsSG
  PUBLIC_HEADERS 
    G4ToolsSGNode.hh
    G4ToolsSGSceneHandler.hh
    G4ToolsSGViewer.hh
  SOURCES
    G4ToolsSGSceneHandler.cc)

geant4_module_link_libraries(G4ToolsSG
  PUBLIC
    G4UIbasic
    G4UIcommon
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

  geant4_module_link_libraries(G4ToolsSG PUBLIC Qt${QT_VERSION_MAJOR}::OpenGL Qt${QT_VERSION_MAJOR}::Gui Qt${QT_VERSION_MAJOR}::PrintSupport Qt${QT_VERSION_MAJOR}::Widgets OpenGL::GL)
endif()

# Windows sources if selected
if(GEANT4_USE_TOOLSSG_WINDOWS_GLES)
  geant4_module_sources(G4ToolsSG PUBLIC_HEADERS G4ToolsSGWindowsGLES.hh SOURCES G4ToolsSGWindowsGLES.cc)
  geant4_module_link_libraries(G4ToolsSG PUBLIC OpenGL::GL)
endif()




