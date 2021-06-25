# - G4gl2ps module build definition

# We need to add definitions depending on what GL drivers are built
if(GEANT4_USE_OPENGL)
  add_definitions(-DG4VIS_BUILD_OPENGL_DRIVER)
endif()

if(GEANT4_USE_INVENTOR)
  add_definitions(-DG4VIS_BUILD_OI_DRIVER)
endif()

# Link to appropriate GL library for the platform/drivers being built
set(G4GL2PS_GL_LIBRARIES OpenGL::GL)
if(APPLE AND (GEANT4_USE_OPENGL_X11 OR GEANT4_USE_INVENTOR_XT OR GEANT4_USE_XM))
  set(G4GL2PS_GL_LIBRARIES XQuartzGL::GL)
endif()

# Define the Geant module
geant4_add_module(G4gl2ps
  PUBLIC_HEADERS
    G4OpenGL2PSAction.hh
    Geant4_gl2ps.h
    gl2ps.h
  SOURCES
    G4OpenGL2PSAction.cc
    gl2ps.cc)

geant4_module_link_libraries(G4gl2ps
  PUBLIC
    G4globman
    ${ZLIB_LIBRARIES}
    ${G4GL2PS_GL_LIBRARIES})

