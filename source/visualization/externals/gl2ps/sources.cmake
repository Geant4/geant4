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
# $Id: sources.cmake 70646 2013-06-03 15:07:58Z gcosmo $
#
#-----------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})
include_directories(${ZLIB_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)

# Must have GL headers available
include_directories(${OPENGL_INCLUDE_DIR})

# We need to add definitions depending on what GL drivers are built
#
if(GEANT4_USE_OPENGL)
  add_definitions(-DG4VIS_BUILD_OPENGL_DRIVER)
endif()

if(GEANT4_USE_INVENTOR)
  add_definitions(-DG4VIS_BUILD_OI_DRIVER)
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
  GLOBAL_DEPENDENCIES
  LINK_LIBRARIES
    ${ZLIB_LIBRARIES}
    ${OPENGL_LIBRARIES}
  )

# List any source specific properties here

