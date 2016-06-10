# - build library

file (GLOB PY_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.py)

foreach (_pyfile ${PY_FILES})
  configure_file (${CMAKE_CURRENT_SOURCE_DIR}/${_pyfile}
                  ${CMAKE_CURRENT_BINARY_DIR}/${_pyfile} COPYONLY)
  list (APPEND PYC_FILES "${CMAKE_CURRENT_BINARY_DIR}/${_pyfile}c")
  list (APPEND PYO_FILES "${CMAKE_CURRENT_BINARY_DIR}/${_pyfile}o")
endforeach()

add_custom_target (Geant4Py ALL)

add_custom_command (
  TARGET Geant4Py
  COMMAND ${PYTHON_EXECUTABLE}
  ARGS -m compileall -f ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND ${PYTHON_EXECUTABLE}
  ARGS -O -m compileall -f ${CMAKE_CURRENT_BINARY_DIR}
)

# install

install (FILES ${PY_FILES} DESTINATION ${G4MODULES_INSTALL_DIR})
install (FILES ${PYC_FILES} DESTINATION ${G4MODULES_INSTALL_DIR})
install (FILES ${PYO_FILES} DESTINATION ${G4MODULES_INSTALL_DIR})
