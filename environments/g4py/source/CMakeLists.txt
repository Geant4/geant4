# - add libs components
# Core Geant4 bindings
add_subdirectory(global)
add_subdirectory(interface)
add_subdirectory(intercoms)
add_subdirectory(run)
add_subdirectory(event)
add_subdirectory(tracking)
add_subdirectory(track)
add_subdirectory(geometry)
add_subdirectory(materials)
add_subdirectory(particles)
add_subdirectory(processes)
add_subdirectory(digits_hits)
add_subdirectory(visualization)
add_subdirectory(graphics_reps)
add_subdirectory(physics_lists)

# Only if GDML
if(GEANT4_USE_GDML)
  add_subdirectory(gdml)
endif()

# Pure Python Geant4 extensions
set(GEANT4_PYTHON_MODULES
  __init__.py
  colortable.py
  g4thread.py
  g4viscp.py
  hepunit.py)
set(GEANT4_PYTHON_BYTECOMPILE_INPUT)
set(GEANT4_PYTHON_BYTECOMPILE_OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/__pycache__")

foreach(_pymod ${GEANT4_PYTHON_MODULES})
  # Generate used to support testing and multimode so .py/pyc files are in the
  # same place as the output extensions
  file(GENERATE 
    OUTPUT ${GEANT4_PYTHON_OUTPUT_DIR}/Geant4/${_pymod}
    INPUT ${CMAKE_CURRENT_SOURCE_DIR}/${_pymod})
  # Configure module for byte-compilation at install time
  # We do this locally because add_custom_command cannot use genexes in
  # its OUTPUT command
  configure_file(${_pymod} ${CMAKE_CURRENT_BINARY_DIR}/${_pymod} COPYONLY)
  list(APPEND GEANT4_PYTHON_BYTECOMPILE_INPUT ${CMAKE_CURRENT_BINARY_DIR}/${_pymod})
endforeach()

add_custom_command(
  OUTPUT "${GEANT4_PYTHON_BYTECOMPILE_OUTPUT}"
  COMMAND ${PYTHON_EXECUTABLE} -m compileall -q ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS ${GEANT4_PYTHON_BYTECOMPILE_INPUT}
  COMMENT "Byte-compiling Geant4 .py modules")

add_custom_target(pyG4ByteCompile ALL DEPENDS "${GEANT4_PYTHON_BYTECOMPILE_OUTPUT}")

install(FILES ${GEANT4_PYTHON_BYTECOMPILE_INPUT} DESTINATION "${CMAKE_INSTALL_PYTHONDIR}/Geant4")
install(DIRECTORY ${GEANT4_PYTHON_BYTECOMPILE_OUTPUT} DESTINATION "${CMAKE_INSTALL_PYTHONDIR}/Geant4")

