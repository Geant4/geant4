#------------------------------------------------------------------------------
# sources.cmake
# Module : G4navigation
# Package: Geant4.src.G4geometry.G4navigation
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 73633 2013-09-03 09:56:54Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/magneticfield/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/volumes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/graphics_reps/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)
include_directories(${CMAKE_SOURCE_DIR}/source/materials/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4navigation
    HEADERS
        G4AuxiliaryNavServices.hh
        G4AuxiliaryNavServices.icc
        G4BrentLocator.hh
        G4DrawVoxels.hh
        G4ErrorPropagationNavigator.hh
        G4GeomTestVolume.hh
        G4GeometryMessenger.hh
        G4GlobalMagFieldMessenger.hh
        G4MultiLevelLocator.hh
        G4MultiNavigator.hh
        G4NavigationLogger.hh
        G4Navigator.hh
        G4Navigator.icc
        G4NormalNavigation.hh
        G4NormalNavigation.icc
        G4ParameterisedNavigation.hh
        G4ParameterisedNavigation.icc
        G4PartialPhantomParameterisation.hh
        G4PathFinder.hh
        G4PhantomParameterisation.hh
        G4PhantomParameterisation.icc
        G4PropagatorInField.hh
        G4PropagatorInField.icc
        G4RegularNavigation.hh
        G4RegularNavigationHelper.hh
        G4ReplicaNavigation.hh
        G4ReplicaNavigation.icc
        G4SafetyHelper.hh
        G4SimpleLocator.hh
        G4TransportationManager.hh
        G4TransportationManager.icc
        G4VIntersectionLocator.hh
        G4VIntersectionLocator.icc
        G4VoxelNavigation.hh
        G4VoxelNavigation.icc
        G4VoxelSafety.hh
    SOURCES
        G4AuxiliaryNavServices.cc
        G4BrentLocator.cc
        G4DrawVoxels.cc
        G4ErrorPropagationNavigator.cc
        G4GeomTestVolume.cc
        G4GeometryMessenger.cc
        G4GlobalMagFieldMessenger.cc
        G4MultiLevelLocator.cc
        G4MultiNavigator.cc
        G4NavigationLogger.cc
        G4Navigator.cc
        G4NormalNavigation.cc
        G4ParameterisedNavigation.cc
        G4PartialPhantomParameterisation.cc
        G4PathFinder.cc
        G4PhantomParameterisation.cc
        G4PropagatorInField.cc
        G4RegularNavigation.cc
        G4RegularNavigationHelper.cc
        G4ReplicaNavigation.cc
        G4SafetyHelper.cc
        G4SimpleLocator.cc
        G4TransportationManager.cc
        G4VIntersectionLocator.cc
        G4VoxelNavigation.cc
        G4VoxelSafety.cc
    GRANULAR_DEPENDENCIES
        G4geometrymng
        G4globman
        G4graphics_reps
        G4intercoms
        G4magneticfield
        G4materials
        G4volumes
    GLOBAL_DEPENDENCIES
        G4global
        G4graphics_reps
        G4intercoms
        G4materials
    LINK_LIBRARIES
)

# List any source specific properties here

