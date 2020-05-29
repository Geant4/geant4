# - build library

# library
geant4_add_pymodule(pyG4run
  pyG4Run.cc
  pyG4RunManager.cc
  pyG4RunManagerKernel.cc
  pyG4UserRunAction.cc
  pyG4VModularPhysicsList.cc
  pyG4VPhysicsConstructor.cc
  pyG4VUserDetectorConstruction.cc
  pyG4VUserPhysicsList.cc
  pyG4VUserPrimaryGeneratorAction.cc
  pymodG4run.cc
)

target_link_libraries(pyG4run PRIVATE G4run G4readout)
install(TARGETS pyG4run DESTINATION "${CMAKE_INSTALL_PYTHONDIR}/Geant4")
