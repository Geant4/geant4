#------------------------------------------------------------------------------
# sources.cmake
# Module : G4GMocren
# Package: Geant4.src.G4visualization.G4GMocren
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 99176 2016-09-07 09:46:36Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})
include_directories(${USOLIDS_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/digits_hits/digits/include)
include_directories(${CMAKE_SOURCE_DIR}/source/digits_hits/hits/include)
include_directories(${CMAKE_SOURCE_DIR}/source/digits_hits/utils/include)
include_directories(${CMAKE_SOURCE_DIR}/source/event/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/navigation/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/CSG/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/specific/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPNumerics/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/graphics_reps/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)
include_directories(${CMAKE_SOURCE_DIR}/source/materials/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/tracking/include)
include_directories(${CMAKE_SOURCE_DIR}/source/visualization/FukuiRenderer/include)
include_directories(${CMAKE_SOURCE_DIR}/source/visualization/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/visualization/modeling/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4GMocren
    HEADERS
        G4GMocrenFile.hh
        G4GMocrenFileCTtoDensityMap.hh
        G4GMocrenFileSceneHandler.hh
        G4GMocrenFileViewer.hh
        G4GMocrenIO.hh
        G4GMocrenMessenger.hh
        G4GMocrenTouchable.hh
    SOURCES
        G4GMocrenFile.cc
        G4GMocrenFileCTtoDensityMap.cc
        G4GMocrenFileSceneHandler.cc
        G4GMocrenFileViewer.cc
        G4GMocrenIO.cc
        G4GMocrenMessenger.cc
    GRANULAR_DEPENDENCIES
        G4FR
        G4csg
        G4detutils
        G4digits
        G4event
        G4geometrymng
        G4globman
        G4graphics_reps
        G4hepnumerics
        G4hits
        G4intercoms
        G4materials
        G4modeling
        G4navigation
        G4partman
        G4specsolids
        G4tracking
        G4vis_management
    GLOBAL_DEPENDENCIES
        G4FR
        G4digits_hits
        G4event
        G4geometry
        G4global
        G4graphics_reps
        G4intercoms
        G4materials
        G4modeling
        G4particles
        G4tracking
        G4vis_management
    LINK_LIBRARIES
)

# List any source specific properties here

