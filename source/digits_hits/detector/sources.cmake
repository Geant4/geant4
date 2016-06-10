#------------------------------------------------------------------------------
# sources.cmake
# Module : G4detector
# Package: Geant4.src.G4digits_hits.G4detector
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 92076 2015-08-17 07:05:45Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/digits_hits/hits/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/biasing/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/navigation/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/volumes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)
include_directories(${CMAKE_SOURCE_DIR}/source/materials/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4detector
    HEADERS
        G4CellScoreComposer.hh
        G4CellScoreValues.hh
        G4CollectionNameVector.hh
        G4HCtable.hh
        G4MultiFunctionalDetector.hh
        G4SDManager.hh
        G4SDStructure.hh
        G4SDmessenger.hh
        G4SensitiveVolumeList.hh
        G4SensitiveVolumeList.icc
        G4TrackLogger.hh
        G4VPrimitiveScorer.hh
        G4VReadOutGeometry.hh
        G4VSDFilter.hh
        G4VSensitiveDetector.hh
	G4MultiSensitiveDetector.hh
    SOURCES
        G4CellScoreComposer.cc
        G4HCtable.cc
        G4MultiFunctionalDetector.cc
        G4SDManager.cc
        G4SDStructure.cc
        G4SDmessenger.cc
        G4SensitiveVolumeList.cc
        G4TrackLogger.cc
        G4VPrimitiveScorer.cc
        G4VReadOutGeometry.cc
        G4VSDFilter.cc
        G4VSensitiveDetector.cc
	G4MultiSensitiveDetector.cc
    GRANULAR_DEPENDENCIES
        G4geombias
        G4geometrymng
        G4globman
        G4hits
        G4intercoms
        G4materials
        G4navigation
        G4partman
        G4track
        G4volumes
    GLOBAL_DEPENDENCIES
        G4geometry
        G4global
        G4intercoms
        G4materials
        G4particles
        G4track
    LINK_LIBRARIES
)

# List any source specific properties here

