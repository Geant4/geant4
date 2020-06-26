#------------------------------------------------------------------------------
# Module : G4volumes
# Package: Geant4.src.G4geometry.G4volumes
#------------------------------------------------------------------------------

#
# Define the Geant4 Module.
#
geant4_define_module(NAME G4volumes
  HEADERS
    G4AssemblyTriplet.hh
    G4AssemblyTriplet.icc
    G4AssemblyStore.hh
    G4AssemblyVolume.hh
    G4AssemblyVolume.icc
    G4EnhancedVecAllocator.hh
    G4GeometryWorkspace.hh
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
	G4VExternalPhysicalVolume.hh
  SOURCES
    G4AssemblyStore.cc
    G4AssemblyVolume.cc
    G4GeometryWorkspace.cc
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
	G4VExternalPhysicalVolume.cc
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
)

# List any source specific properties here
