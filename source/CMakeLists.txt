# - Define toolkit source code modules and composed libraries

add_subdirectory(analysis)
add_subdirectory(digits_hits)
add_subdirectory(error_propagation)
add_subdirectory(event)
add_subdirectory(externals)

# Target G4g3tog4 is optional based on user selection...
if(GEANT4_USE_G3TOG4)
  add_subdirectory(g3tog4)
endif()

add_subdirectory(geometry)
add_subdirectory(global)
add_subdirectory(graphics_reps)
add_subdirectory(intercoms)
add_subdirectory(interfaces)
add_subdirectory(materials)
add_subdirectory(parameterisations)
add_subdirectory(particles)
add_subdirectory(persistency)
add_subdirectory(physics_lists)
add_subdirectory(processes)
add_subdirectory(readout)
add_subdirectory(run)
add_subdirectory(track)
add_subdirectory(tracking)
add_subdirectory(visualization)

# Optional reset of categories
# This is *only for build testing purposes* as the resultant libraries are not
# usable by tests, examples, or any other client application
if(__GEANT4_LIBRARY_DEFINITION_FILE AND (EXISTS ${PROJECT_SOURCE_DIR}/tests))
  get_filename_component(__lib_def_file "${__GEANT4_LIBRARY_DEFINITION_FILE}" ABSOLUTE)
  if(NOT EXISTS "${__lib_def_file}")
    message(FATAL_ERROR "library definition file '${__lib_def_file}' does not exist")
  endif()

  __geant4_category_reset()
  include("${__lib_def_file}")
endif()

# - Compose libs
geant4_compose_targets()
