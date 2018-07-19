#------------------------------------------------------------------------------
# sources.cmake
# Module : G4detutils
# Package: Geant4.src.G4digits_hits.G4detutils
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 99154 2016-09-07 08:06:30Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/digits_hits/detector/include)
include_directories(${CMAKE_SOURCE_DIR}/source/digits_hits/digits/include)
include_directories(${CMAKE_SOURCE_DIR}/source/digits_hits/hits/include)
include_directories(${CMAKE_SOURCE_DIR}/source/digits_hits/scorer/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/divisions/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/navigation/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/CSG/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/volumes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPNumerics/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/graphics_reps/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)
include_directories(${CMAKE_SOURCE_DIR}/source/materials/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4detutils
    HEADERS
        G4DefaultLinearColorMap.hh
        G4ScoreLogColorMap.hh
        G4ScoreQuantityMessenger.hh
        G4ScoringBox.hh
        G4ScoringCylinder.hh
        G4ScoringManager.hh
        G4ScoringMessenger.hh
        G4VScoreColorMap.hh
        G4VScoreWriter.hh
        G4VScoringMesh.hh
    SOURCES
        G4DefaultLinearColorMap.cc
        G4ScoreLogColorMap.cc
        G4ScoreQuantityMessenger.cc
        G4ScoreQuantityMessengerQCmd.cc
        G4ScoringBox.cc
        G4ScoringCylinder.cc
        G4ScoringManager.cc
        G4ScoringMessenger.cc
        G4VScoreColorMap.cc
        G4VScoreWriter.cc
        G4VScoringMesh.cc
    GRANULAR_DEPENDENCIES
        G4csg
        G4detector
        G4detscorer
        G4digits
        G4geomdivision
        G4geometrymng
        G4globman
        G4graphics_reps
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
        G4graphics_reps
        G4intercoms
        G4materials
        G4particles
        G4track
    LINK_LIBRARIES
)

# List any source specific properties here

