# - build library
set(_TARGET _ezgeom)
g4pytest_add_module(${_TARGET}
  EzDetectorConstruction.cc
  G4EzVolume.cc
  G4EzVoxelParameterization.cc
  G4EzWorld.cc
  pyEzgeom.cc
)
target_link_libraries(${_TARGET} PRIVATE G4materials G4geometry G4run)
