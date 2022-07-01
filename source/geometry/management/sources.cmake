# - G4geometrymng module build definition

# Configure header for preprocessor symbols for USolids
configure_file(${CMAKE_CURRENT_LIST_DIR}/include/G4GeomConfig.hh.in
  ${CMAKE_CURRENT_BINARY_DIR}/include/G4GeomConfig.hh)

# Define the Geant4 Module.
geant4_add_module(G4geometrymng
  PUBLIC_HEADERS
    ${CMAKE_CURRENT_BINARY_DIR}/include/G4GeomConfig.hh
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
    G4GeomTypes.hh
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
    G4VCurvedTrajectoryFilter.cc
    G4VNestedParameterisation.cc
    G4VPVDivisionFactory.cc
    G4VPVParameterisation.cc
    G4VPhysicalVolume.cc
    G4VSolid.cc
    G4VTouchable.cc
    G4VoxelLimits.cc)

# - Add path to generated header
geant4_module_include_directories(G4geometrymng
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>)

# - Link to modules/deps
geant4_module_link_libraries(G4geometrymng
  PUBLIC G4globman G4hepgeometry G4graphics_reps ${VECGEOM_LIBRARIES}
  PRIVATE G4materials G4heprandom)

# WORKAROUND: When building/testing examples uing ROOT, ROOT's
# dictionary generation is not smart enough to handle target usage
# requirements for include paths. Explicitly add the path to the
# generated header into build time include paths...
set_property(GLOBAL APPEND
  PROPERTY GEANT4_BUILDTREE_INCLUDE_DIRS "${CMAKE_CURRENT_BINARY_DIR}/include")

