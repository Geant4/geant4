#------------------------------------------------------------------------------
# sources.cmake
# Module : G4Tree
# Package: Geant4.src.G4visualization.G4Tree
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 88190 2015-02-02 17:24:54Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})
include_directories(${USOLIDS_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/digits_hits/detector/include)
include_directories(${CMAKE_SOURCE_DIR}/source/digits_hits/hits/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/navigation/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/CSG/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/specific/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/volumes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/graphics_reps/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)
include_directories(${CMAKE_SOURCE_DIR}/source/materials/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)
include_directories(${CMAKE_SOURCE_DIR}/source/visualization/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/visualization/modeling/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4Tree
    HEADERS
        G4ASCIITree.hh
        G4ASCIITreeMessenger.hh
        G4ASCIITreeSceneHandler.hh
        G4ASCIITreeViewer.hh
        G4VTree.hh
        G4VTreeSceneHandler.hh
        G4VTreeSceneHandler.icc
        G4VTreeViewer.hh
    SOURCES
        G4ASCIITree.cc
        G4ASCIITreeMessenger.cc
        G4ASCIITreeSceneHandler.cc
        G4ASCIITreeViewer.cc
        G4VTree.cc
        G4VTreeSceneHandler.cc
        G4VTreeViewer.cc
    GRANULAR_DEPENDENCIES
        G4csg
        G4detector
        G4geometrymng
        G4globman
        G4graphics_reps
        G4hits
        G4intercoms
        G4materials
        G4modeling
        G4navigation
        G4partman
        G4specsolids
        G4track
        G4vis_management
        G4volumes
    GLOBAL_DEPENDENCIES
        G4digits_hits
        G4geometry
        G4global
        G4graphics_reps
        G4intercoms
        G4materials
        G4modeling
        G4particles
        G4track
        G4vis_management
    LINK_LIBRARIES
)

# List any source specific properties here

