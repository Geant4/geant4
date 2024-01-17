# - G4volumes module build definition

#
# Define the Geant4 Module.
#
geant4_add_module(G4volumes
  PUBLIC_HEADERS
    G4AssemblyTriplet.hh
    G4AssemblyTriplet.icc
    G4AssemblyStore.hh
    G4AssemblyVolume.hh
    G4AssemblyVolume.icc
    G4EnhancedVecAllocator.hh
    G4GeometryWorkspace.hh
    G4LogicalBorderSurface.hh
    G4LogicalBorderSurface.icc
    G4LogicalSkinSurface.hh
    G4LogicalSkinSurface.icc
    G4PVParameterised.hh
    G4PVPlacement.hh
    G4PVReplica.hh
    G4ReflectionFactory.hh
    G4VExternalPhysicalVolume.hh
  SOURCES
    G4AssemblyStore.cc
    G4AssemblyVolume.cc
    G4GeometryWorkspace.cc
    G4LogicalBorderSurface.cc
    G4LogicalSkinSurface.cc
    G4PVParameterised.cc
    G4PVPlacement.cc
    G4PVReplica.cc
    G4ReflectionFactory.cc
    G4VExternalPhysicalVolume.cc)

geant4_module_link_libraries(G4volumes
  PUBLIC G4globman G4hepgeometry G4geometrymng)
