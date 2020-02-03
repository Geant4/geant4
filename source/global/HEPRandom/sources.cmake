#------------------------------------------------------------------------------
# sources.cmake
# Module : 
# Package: Geant4.src.G4global.
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
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)

#
# Define the Geant4 Module.
#

include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4heprandom 
    HEADERS
        G4Poisson.hh
        G4QuickRand.hh
        G4RandomDirection.hh
        G4RandomTools.hh
	G4UniformRandPool.hh
        Randomize.hh
    SOURCES
        G4Poisson.cc
	G4UniformRandPool.cc
    GRANULAR_DEPENDENCIES
    GLOBAL_DEPENDENCIES
    LINK_LIBRARIES
)

# List any source specific properties here

