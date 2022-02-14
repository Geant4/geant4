# - G4ToolsSG module build definition

# Module has optional sources
# List those always built
set(G4VIS_MODULE_TOOLSSG_HEADERS
  G4ToolsSGNode.hh
  G4ToolsSGSceneHandler.hh
  G4ToolsSGViewer.hh
)

set(G4VIS_MODULE_TOOLSSG_SOURCES
  G4ToolsSGSceneHandler.cc
)

# X11 sources/links if selected
if(GEANT4_USE_TOOLSSG_X11_GLES)
  list(APPEND G4VIS_MODULE_TOOLSSG_HEADERS
    G4ToolsSGX11GLES.hh
  )

  list(APPEND G4VIS_MODULE_TOOLSSG_SOURCES
    G4ToolsSGX11GLES.cc
  )

  list(APPEND G4VIS_MODULE_TOOLSSG_LINK_LIBRARIES X11::SM X11::ICE X11::X11 X11::Xext X11::Xmu)
  if(APPLE)
    list(APPEND G4VIS_MODULE_TOOLSSG_LINK_LIBRARIES XQuartzGL::GL)
  else()
    list(APPEND G4VIS_MODULE_TOOLSSG_LINK_LIBRARIES OpenGL::GL)
  endif()
endif()

# Xt sources/links if selected
if(GEANT4_USE_TOOLSSG_XT_GLES)
  list(APPEND G4VIS_MODULE_TOOLSSG_HEADERS
    G4ToolsSGXtGLES.hh
  )

  list(APPEND G4VIS_MODULE_TOOLSSG_SOURCES
    G4ToolsSGXtGLES.cc
  )

  add_definitions(-DTOOLS_USE_GL_GL_H)

  list(APPEND G4VIS_MODULE_TOOLSSG_LINK_LIBRARIES X11::SM X11::ICE X11::X11 X11::Xext X11::Xmu ${G4VIS_MODULE_TOOLSSG_LINK_LIBRARIES})
  if(APPLE)
    list(APPEND G4VIS_MODULE_TOOLSSG_LINK_LIBRARIES XQuartzGL::GL)
  else()
    list(APPEND G4VIS_MODULE_TOOLSSG_LINK_LIBRARIES OpenGL::GL)
  endif()
endif()

# Qt sources/links if selected
if(GEANT4_USE_TOOLSSG_QT_GLES)
  list(APPEND G4VIS_MODULE_TOOLSSG_HEADERS
    G4ToolsSGQtViewer.hh
    G4ToolsSGQtGLES.hh
  )

  list(APPEND G4VIS_MODULE_TOOLSSG_SOURCES
    G4ToolsSGQtGLES.cc
    G4ToolsSGQtViewer.cc
  )

  list(APPEND G4VIS_MODULE_TOOLSSG_LINK_LIBRARIES Qt5::OpenGL Qt5::Gui Qt5::PrintSupport Qt5::Widgets OpenGL::GL)
endif()


# Windows sources if selected
if(GEANT4_USE_TOOLSSG_WINDOWS_GLES)
  list(APPEND G4VIS_MODULE_TOOLSSG_HEADERS
    G4ToolsSGWindowsGLES.hh
  )

  list(APPEND G4VIS_MODULE_TOOLSSG_SOURCES
    G4ToolsSGWindowsGLES.cc
  )

  list(APPEND G4VIS_MODULE_TOOLSSG_LINK_LIBRARIES OpenGL::GL)
endif()


# Freetype support if selected
if(GEANT4_USE_FREETYPE)
  add_definitions(-DTOOLS_USE_FREETYPE)
  list(APPEND G4VIS_MODULE_TOOLSSG_LINK_LIBRARIES Freetype::Freetype)
endif()

list(REMOVE_DUPLICATES G4VIS_MODULE_TOOLSSG_LINK_LIBRARIES)

# Define the Geant4 Module.
geant4_add_module(G4ToolsSG
  PUBLIC_HEADERS ${G4VIS_MODULE_TOOLSSG_HEADERS}
  SOURCES ${G4VIS_MODULE_TOOLSSG_SOURCES})

geant4_module_link_libraries(G4ToolsSG
  PUBLIC
    G4intercoms
    G4interfaces
    G4modeling
    G4vis_management
    G4tools
    ${G4VIS_MODULE_TOOLSSG_LINK_LIBRARIES}
  PRIVATE
    G4navigation)
