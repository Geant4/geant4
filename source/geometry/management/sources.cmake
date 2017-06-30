#------------------------------------------------------------------------------
# sources.cmake
# Module : G4geometrymng
# Package: Geant4.src.G4geometry.G4geometrymng
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 104304 2017-05-24 08:57:31Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})
include_directories(${USOLIDS_INCLUDE_DIRS})

# List internal includes needed.
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
GEANT4_DEFINE_MODULE(NAME G4geometrymng
    HEADERS
        G4AffineTransform.hh
        G4AffineTransform.icc
        G4BlockingList.hh
        G4BlockingList.icc
        G4BoundingEnvelope.hh
        G4ErrorCylSurfaceTarget.hh
        G4ErrorPlaneSurfaceTarget.hh
        G4ErrorSurfaceTarget.hh
        G4ErrorTanPlaneTarget.hh
        G4ErrorTarget.hh
        G4GeomSplitter.hh
        G4GeomTools.hh
        G4GeometryManager.hh
        G4IdentityTrajectoryFilter.hh
        G4LogicalCrystalVolume.hh
        G4LogicalSurface.hh
        G4LogicalSurface.icc
        G4LogicalVolume.hh
        G4LogicalVolume.icc
        G4LogicalVolumeStore.hh
        G4PhysicalVolumeStore.hh
        G4ReflectedSolid.hh
        G4Region.hh
        G4Region.icc
        G4RegionStore.hh
        G4ScaleTransform.hh
        G4ScaleTransform.icc
        G4SmartVoxelHeader.hh
        G4SmartVoxelHeader.icc
        G4SmartVoxelNode.hh
        G4SmartVoxelNode.icc
        G4SmartVoxelProxy.hh
        G4SmartVoxelProxy.icc
        G4SmartVoxelStat.hh
        G4SolidStore.hh
        G4TouchableHandle.hh
        G4UAdapter.hh
        G4USolid.hh
        G4VCurvedTrajectoryFilter.hh
        G4VNestedParameterisation.hh
        G4VPVDivisionFactory.hh
        G4VPVParameterisation.hh
        G4VPhysicalVolume.hh
        G4VPhysicalVolume.icc
        G4VSolid.hh
        G4VSolid.icc
        G4VStoreNotifier.hh
        G4VTouchable.hh
        G4VTouchable.icc
        G4VUserRegionInformation.hh
        G4VVolumeMaterialScanner.hh
        G4VoxelLimits.hh
        G4VoxelLimits.icc
        geomwdefs.hh
        meshdefs.hh
        voxeldefs.hh
    SOURCES
        G4BlockingList.cc
        G4BoundingEnvelope.cc
        G4ErrorCylSurfaceTarget.cc
        G4ErrorPlaneSurfaceTarget.cc
        G4ErrorSurfaceTarget.cc
        G4ErrorTanPlaneTarget.cc
        G4ErrorTarget.cc
        G4GeomTools.cc
        G4GeometryManager.cc
        G4IdentityTrajectoryFilter.cc
        G4LogicalCrystalVolume.cc
        G4LogicalSurface.cc
        G4LogicalVolume.cc
        G4LogicalVolumeStore.cc
        G4PhysicalVolumeStore.cc
        G4ReflectedSolid.cc
        G4Region.cc
        G4RegionStore.cc
        G4SmartVoxelHeader.cc
        G4SmartVoxelNode.cc
        G4SmartVoxelProxy.cc
        G4SmartVoxelStat.cc
        G4SolidStore.cc
        G4USolid.cc
        G4VCurvedTrajectoryFilter.cc
        G4VNestedParameterisation.cc
        G4VPVDivisionFactory.cc
        G4VPVParameterisation.cc
        G4VPhysicalVolume.cc
        G4VSolid.cc
        G4VTouchable.cc
        G4VoxelLimits.cc
    GRANULAR_DEPENDENCIES
        G4globman
        G4graphics_reps
        G4intercoms
        G4materials
    GLOBAL_DEPENDENCIES
        G4global
        G4graphics_reps
        G4intercoms
        G4materials
    LINK_LIBRARIES
        ${USOLIDS_LIBRARIES}
)

# List any source specific properties here

