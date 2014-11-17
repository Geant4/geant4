# - build library

# library
set(_TARGET pyG4particles)
add_library(
  ${_TARGET} SHARED
  pyG4DecayTable.cc
  pyG4DynamicParticle.cc
  pyG4ParticleDefinition.cc
  pyG4ParticleList.cc
  pyG4ParticleTable.cc
  pyG4PrimaryParticle.cc
  pyG4PrimaryVertex.cc
  pymodG4particles.cc
)

set_target_properties(${_TARGET} PROPERTIES PREFIX "")
set_target_properties(${_TARGET} PROPERTIES OUTPUT_NAME "G4particles")
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
