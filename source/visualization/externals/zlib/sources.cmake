#------------------------------------------------------------------------------
# sources.cmake
# Module : G4zlib
# Package: Geant4.src.G4visualization..G4zlib
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake,v 1.1 2010-09-29 19:14:24 bmorgan Exp $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4zlib
    HEADERS
        crc32.h
        deflate.h
        trees.h
        zconf.h
        zlib.h
        zutil.h
    SOURCES
        adler32.cc
        compress.cc
        crc32.cc
        deflate.cc
        trees.cc
        zutil.cc
    GRANULAR_DEPENDENCIES
    GLOBAL_DEPENDENCIES
    LINK_LIBRARIES
)

# List any source specific properties here

