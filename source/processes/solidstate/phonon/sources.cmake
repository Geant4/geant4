#------------------------------------------------------------------------------
# sources.cmake
# Module : G4phonon
# Package: Geant4.src.G4processes.G4phonon
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/10/2013
#
# $Id: sources.cmake 75725 2013-11-05 16:52:30Z mkelsey $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/volumes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/materials/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/bosons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4phonon
    HEADERS
	G4LatticeManager.hh
	G4LatticeReader.hh
	G4PhononDownconversion.hh
	G4PhononPolarization.hh
	G4PhononReflection.hh
	G4PhononScattering.hh
	G4PhononTrackMap.hh
	G4VPhononProcess.hh
    SOURCES
	G4LatticeManager.cc
	G4LatticeReader.cc
	G4PhononDownconversion.cc
	G4PhononPolarization.cc
	G4PhononReflection.cc
	G4PhononScattering.cc
	G4PhononTrackMap.cc
	G4VPhononProcess.cc
    GRANULAR_DEPENDENCIES
        G4bosons
        G4geometrymng
        G4globman
        G4materials
        G4partman
        G4procman
        G4track
        G4volumes
    GLOBAL_DEPENDENCIES
        G4geometry
        G4global
        G4materials
        G4particles
        G4track
    LINK_LIBRARIES
)

# List any source specific properties here

