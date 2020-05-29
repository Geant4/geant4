# - build test library

# library
set(_TARGET gtest01)
g4py_add_module(${_TARGET}
  Particles.cc
  PhysicsListEMstd.cc
  QDetectorConstruction.cc
  QEventAction.cc
  QMaterials.cc
  QPhysicsList.cc
  QPrimaryGeneratorAction.cc
  gtest01.cc
)
target_link_libraries(${_TARGET} PRIVATE G4particles G4processes G4run)

