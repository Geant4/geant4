#-----------------------------------------------------------------------
# Module : G4gl2ps
# Package: Geant4.src.G4visualization..G4gl2ps
#-----------------------------------------------------------------------

# We need to add definitions depending on what GL drivers are built
#
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

# Define the Geant4 Module.
#
GEANT4_DEFINE_MODULE(NAME G4gl2ps
  HEADERS
    G4OpenGL2PSAction.hh
    Geant4_gl2ps.h
    gl2ps.h
  SOURCES
    G4OpenGL2PSAction.cc
    gl2ps.cc
  GRANULAR_DEPENDENCIES
    G4globman
  GLOBAL_DEPENDENCIES
    G4global
  LINK_LIBRARIES
    ${ZLIB_LIBRARIES}
    ${G4GL2PS_GL_LIBRARIES}
  SOURCES_EXCLUDE_FORMAT
    gl2ps.h
  )

# List any source specific properties here

