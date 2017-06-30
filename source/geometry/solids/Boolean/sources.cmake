#------------------------------------------------------------------------------
# sources.cmake
# Module : G4geomBoolean
# Package: Geant4.src.G4geometry.G4geomBoolean
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 103464 2017-04-11 07:24:03Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})
include_directories(${USOLIDS_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/volumes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/CSG/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/specific/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/graphics_reps/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4geomBoolean
    HEADERS
        G4BooleanSolid.hh
        G4BooleanSolid.icc
        G4DisplacedSolid.hh
        G4IntersectionSolid.hh
        G4MultiUnion.hh
        G4ScaledSolid.hh
        G4SubtractionSolid.hh
        G4UnionSolid.hh
    SOURCES
        G4BooleanSolid.cc
        G4DisplacedSolid.cc
        G4IntersectionSolid.cc
        G4MultiUnion.cc
        G4ScaledSolid.cc
        G4SubtractionSolid.cc
        G4UnionSolid.cc
    GRANULAR_DEPENDENCIES
        G4geometrymng
        G4globman
        G4graphics_reps
        G4intercoms
        G4volumes
        G4csg
        G4specsolids
    GLOBAL_DEPENDENCIES
        G4global
        G4graphics_reps
        G4intercoms
    LINK_LIBRARIES
        ${USOLIDS_LIBRARIES}
)

# List any source specific properties here

