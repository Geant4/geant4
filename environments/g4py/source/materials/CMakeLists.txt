# - build library

# library
set(_TARGET pyG4materials)
add_library(
  ${_TARGET} SHARED
  pyG4AtomicShells.cc
  pyG4Element.cc
  pyG4ElementTable.cc
  pyG4ElementVector.cc
  pyG4Isotope.cc
  pyG4Material.cc
  pyG4MaterialTable.cc
  pyG4NistManager.cc
  pymodG4materials.cc
)

set_target_properties(${_TARGET} PROPERTIES PREFIX "")
set_target_properties(${_TARGET} PROPERTIES OUTPUT_NAME "G4materials")
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
