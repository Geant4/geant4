#------------------------------------------------------------------------------
# sources.cmake
# Module : G4volumes
# Package: Geant4.src.G4geometry.G4volumes
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 80876 2014-05-14 08:51:07Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/magneticfield/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/graphics_reps/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4volumes
    HEADERS
        G4AssemblyTriplet.hh
        G4AssemblyTriplet.icc
        G4AssemblyVolume.hh
        G4AssemblyVolume.icc
        G4EnhancedVecAllocator.hh
        G4GeometryWorkspace.hh
        G4GeometryWorkspacePool.hh
        G4GRSSolid.hh
        G4GRSSolid.icc
        G4GRSSolidHandle.hh
        G4GRSVolume.hh
        G4GRSVolume.icc
        G4GRSVolumeHandle.hh
        G4LogicalBorderSurface.hh
        G4LogicalBorderSurface.icc
        G4LogicalSkinSurface.hh
        G4LogicalSkinSurface.icc
        G4NavigationHistory.hh
        G4NavigationHistory.icc
        G4NavigationHistoryPool.hh
        G4NavigationLevel.hh
        G4NavigationLevel.icc
        G4NavigationLevelRep.hh
        G4NavigationLevelRep.icc
        G4PVParameterised.hh
        G4PVPlacement.hh
        G4PVReplica.hh
        G4ReflectionFactory.hh
        G4TouchableHistory.hh
        G4TouchableHistory.icc
        G4TouchableHistoryHandle.hh
    SOURCES
        G4AssemblyVolume.cc
        G4GeometryWorkspace.cc
        G4GeometryWorkspacePool.cc
        G4GRSSolid.cc
        G4GRSVolume.cc
        G4LogicalBorderSurface.cc
        G4LogicalSkinSurface.cc
        G4NavigationHistory.cc
        G4NavigationHistoryPool.cc
        G4NavigationLevel.cc
        G4NavigationLevelRep.cc
        G4PVParameterised.cc
        G4PVPlacement.cc
        G4PVReplica.cc
        G4ReflectionFactory.cc
        G4TouchableHistory.cc
    GRANULAR_DEPENDENCIES
        G4geometrymng
        G4globman
        G4graphics_reps
        G4intercoms
        G4magneticfield
    GLOBAL_DEPENDENCIES
        G4global
        G4graphics_reps
        G4intercoms
    LINK_LIBRARIES
)

# List any source specific properties here

