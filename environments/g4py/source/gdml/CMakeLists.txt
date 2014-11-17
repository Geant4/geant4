# - build library

if(GEANT4_HAS_GDML)
  add_definitions(-DENABLE_GDML)
endif()

include_directories (${XERCESC_INCLUDE_DIRS})

# library
set(_TARGET pyG4gdml)
add_library(
  ${_TARGET} SHARED
  pyG4GDMLParser.cc
  pymodG4g4gdml.cc
)

set_target_properties(${_TARGET} PROPERTIES PREFIX "")
set_target_properties(${_TARGET} PROPERTIES OUTPUT_NAME "G4gdml")
set_target_properties(${_TARGET} PROPERTIES SUFFIX ".so")
set_target_properties(${_TARGET}
                      PROPERTIES INSTALL_RPATH
                      ${GEANT4_LIBRARY_DIR}
                      BUILD_WITH_INSTALL_RPATH TRUE)

target_link_libraries (${_TARGET}
                      ${GEANT4_LIBRARIES_WITH_VIS} ${XERCESC_LIBRARY}
                      boost_python ${PYTHON_LIBRARIES})

# install
install(TARGETS ${_TARGET} LIBRARY DESTINATION ${G4MODULES_INSTALL_DIR})
