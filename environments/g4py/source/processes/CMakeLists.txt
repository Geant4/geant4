# - build library

# library
set(_TARGET pyG4processes)
add_library(
  ${_TARGET} SHARED
  pyG4CrossSectionHandler.cc
  pyG4EmCalculator.cc
  pyG4LossTableManager.cc
  pyG4ProcVector.cc
  pyG4ProcessManager.cc
  pyG4ProcessTable.cc
  pyG4ProcessType.cc
  pyG4ProductionCutsTable.cc
  pyG4VCrossSectionHandler.cc
  pyG4VProcess.cc
  pymodG4processes.cc
)

set_target_properties(${_TARGET} PROPERTIES PREFIX "")
set_target_properties(${_TARGET} PROPERTIES OUTPUT_NAME "G4processes")
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
