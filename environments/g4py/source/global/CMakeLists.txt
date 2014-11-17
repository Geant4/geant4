# - build library

# library
set(_TARGET pyG4global)
add_library(
  ${_TARGET} SHARED
  G4PyCoutDestination.cc
  pyG4ApplicationState.cc
  pyG4Exception.cc
  pyG4ExceptionHandler.cc
  pyG4ExceptionSeverity.cc
  pyG4RandomDirection.cc
  pyG4RotationMatrix.cc
  pyG4StateManager.cc
  pyG4String.cc
  pyG4ThreeVector.cc
  pyG4Timer.cc
  pyG4Transform3D.cc
  pyG4TwoVector.cc
  pyG4UnitsTable.cc
  pyG4UserLimits.cc
  pyG4Version.cc
  pyRandomEngines.cc
  pyRandomize.cc
  pygeomdefs.cc
  pyglobals.cc
  pymodG4global.cc
)

set_target_properties(${_TARGET} PROPERTIES PREFIX "")
set_target_properties(${_TARGET} PROPERTIES OUTPUT_NAME "G4global")
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

