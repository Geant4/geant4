# - build library

# library
set(_TARGET pyG4run)
add_library(
  ${_TARGET} SHARED
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

set_target_properties(${_TARGET} PROPERTIES PREFIX "")
set_target_properties(${_TARGET} PROPERTIES OUTPUT_NAME "G4run")
set_target_properties(${_TARGET} PROPERTIES SUFFIX ".so")
set_target_properties(${_TARGET}
                      PROPERTIES INSTALL_RPATH
                      ${GEANT4_LIBRARY_DIR}
                      BUILD_WITH_INSTALL_RPATH TRUE)

target_link_libraries (${_TARGET}
                      ${GEANT4_LIBRARIES_WITH_VIS} boost_python
                      ${PYTHON_LIBRARIES})

# install
install(TARGETS ${_TARGET} LIBRARY DESTINATION ${G4MODULES_INSTALL_DIR})
