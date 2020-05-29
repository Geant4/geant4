# Try and get modules output to config dependent site-packages
# Can put all Python outputs under GEANT4_PYTHON_OUTPUT_DIR (and conf variants)
# That's the baseline, packages however need to then adjust output path for
# any substructure
set(GEANT4_PYTHON_OUTPUT_DIR "$<TARGET_FILE_DIR:G4global>/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages")

set(CMAKE_INSTALL_PYTHONDIR "${CMAKE_INSTALL_LIBDIR}/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages")
if(NOT IS_ABSOLUTE "${CMAKE_INSTALL_PYTHONDIR}")
  set(CMAKE_INSTALL_FULL_PYTHONDIR "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_PYTHONDIR}")
else()
  set(CMAKE_INSTALL_FULL_PYTHONDIR "${CMAKE_INSTALL_PYTHONDIR}")
endif()

# Extension modules link to libraries so relative RPATH modules to libraries
file(RELATIVE_PATH GEANT4_PYMODULE_RPATH
  "${CMAKE_INSTALL_FULL_PYTHONDIR}/Geant4"
  "${CMAKE_INSTALL_FULL_LIBDIR}")

# Basic "add module" command to wrap common functionality
# Should also help when migrating to pybind11 or FindPython builtins
function(g4py_add_module target_name)
  add_library(${target_name} MODULE ${ARGN})

  # Python extension module naming
  set_property(TARGET ${target_name} PROPERTY PREFIX "")
  if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    set_property(TARGET ${target_name} PROPERTY SUFFIX ".pyd")
    set_property(TARGET ${target_name} PROPERTY DEBUG_SUFFIX "_d")
  endif()

  # Adjust output location (may need adjustment on Windows)
  set_property(TARGET ${target_name} PROPERTY LIBRARY_OUTPUT_DIRECTORY "${GEANT4_PYTHON_OUTPUT_DIR}")
  foreach(_conftype ${CMAKE_CONFIGURATION_TYPES})
    string(TOUPPER ${_conftype} _conftype_uppercase)
    set_property(TARGET ${target_name} PROPERTY LIBRARY_OUTPUT_DIRECTORY_${_conftype_uppercase} "${GEANT4_PYTHON_OUTPUT_DIR}")
  endforeach()

  # Force linkage to extension system, other deps to be linked by client
  # Use dynamic_lookup on Darwin to avoid direct linking to libpython (Linux has this as default,
  # to be checked for Windows)
  target_include_directories(${target_name} PRIVATE ${PYTHON_INCLUDE_DIRS})
  target_link_libraries(${target_name} PRIVATE Boost::python "$<$<PLATFORM_ID:Darwin>:-undefined dynamic_lookup>")
endfunction()


# Handle core Geant4 modules
# Wraps g4py_add_module, because we have to avoid duplicate target names
# (e.g. "G4global" for libG4global.dylib vs G4global.so)
# It just manages the naming and output name/directory
function(geant4_add_pymodule _name)
  if(NOT (_name MATCHES "^pyG4[a-zA-Z0-9_]+"))
    message(FATAL_ERROR "Invalid Geant4 Python module name '${_name}'. Names must begin with 'pyG4'")
  endif()

  set(_module_target "${_name}")
  string(REGEX REPLACE "^py" "" _module_lib "${_module_target}")

  g4py_add_module(${_module_target} ${ARGN})

  set_property(TARGET ${_module_target} PROPERTY OUTPUT_NAME "${_module_lib}")

  set_property(TARGET ${_module_target} APPEND_STRING PROPERTY LIBRARY_OUTPUT_DIRECTORY "/Geant4")
  foreach(_conftype ${CMAKE_CONFIGURATION_TYPES})
    string(TOUPPER ${_conftype} _conftype_uppercase)
    set_property(TARGET ${_module_target} APPEND_STRING PROPERTY LIBRARY_OUTPUT_DIRECTORY_${_conftype_uppercase} "/Geant4")
  endforeach()

  if(UNIX AND NOT APPLE)
    set_property(TARGET ${_module_target} PROPERTY
      INSTALL_RPATH "\$ORIGIN/${GEANT4_PYMODULE_RPATH}")
  elseif(APPLE)
    set_property(TARGET ${_module_target} PROPERTY
      INSTALL_RPATH "@loader_path/${GEANT4_PYMODULE_RPATH}")
  endif()
endfunction()
