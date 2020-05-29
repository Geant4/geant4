# - Function to build the C++ modules
function(g4pytest_add_module _TARGET)
  g4py_add_module(${_TARGET} ${ARGN})

  # Filthy temp hack for submodule paths....
  string(REGEX REPLACE "^_" "" _SUBMODULE "${_TARGET}")
  set_property(TARGET ${_TARGET} APPEND_STRING PROPERTY LIBRARY_OUTPUT_DIRECTORY "/g4pytest/${_SUBMODULE}")
  foreach(_conftype ${CMAKE_CONFIGURATION_TYPES})
    string(TOUPPER ${_conftype} _conftype_uppercase)
    set_property(TARGET ${_TARGET} APPEND_STRING PROPERTY LIBRARY_OUTPUT_DIRECTORY_${_conftype_uppercase} "/g4pytest/${_SUBMODULE}")
  endforeach()
endfunction()


# - add libs components
add_subdirectory(ExN01geom)
add_subdirectory(ExN03geom)
add_subdirectory(Qgeom)
add_subdirectory(ezgeom)
add_subdirectory(NISTmaterials)
add_subdirectory(Qmaterials)
add_subdirectory(EMSTDpl)
add_subdirectory(ExN01pl)
add_subdirectory(MedicalBeam)
add_subdirectory(ParticleGun)
# Remainder are pure python

# Copy/configure pure python components
file(GLOB_RECURSE PY_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.py)
foreach(_pyfile ${PY_FILES})
  file(GENERATE 
    OUTPUT ${GEANT4_PYTHON_OUTPUT_DIR}/g4pytest/${_pyfile}
    INPUT ${CMAKE_CURRENT_SOURCE_DIR}/${_pyfile}
    )
endforeach()
