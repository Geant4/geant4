#-----------------------------------------------------------------------
# sources.cmake
# Module : G4gl2ps
# Package: Geant4.src.G4visualization..G4gl2ps
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
#
#-----------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})
include_directories(${ZLIB_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${PROJECT_SOURCE_DIR}/source/global/management/include)
# WORKAROUND: Now have a generated header, so must add include directory,
# but not that we lose one directory level because the management subcategory
# is merged into the main global one!
include_directories(${PROJECT_BINARY_DIR}/source/global/include)

# We need to add definitions depending on what GL drivers are built
#
if(GEANT4_USE_OPENGL)
  add_definitions(-DG4VIS_BUILD_OPENGL_DRIVER)
endif()

if(GEANT4_USE_INVENTOR)
  add_definitions(-DG4VIS_BUILD_OI_DRIVER)
endif()

set(G4GL2PS_GL_LIBRARIES OpenGL::GL OpenGL::GLU)
if(APPLE AND GEANT4_USE_OPENGL_X11)
  list(APPEND G4GL2PS_GL_LIBRARIES XQuartzGL::GL XQuartzGL::GLU)
endif()

# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
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

