#------------------------------------------------------------------------------
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
# $Id: sources.cmake,v 1.1 2010-09-29 19:14:08 bmorgan Exp $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/visualization/externals/zlib/include)

#
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
        G4zlib
    GLOBAL_DEPENDENCIES
    LINK_LIBRARIES
)

# List any source specific properties here

