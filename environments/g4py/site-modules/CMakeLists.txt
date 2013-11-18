# - add libs components

set (G4SITEMODULES_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR}/g4py )

include_directories (
  ${PYTHON_INCLUDE_PATH}
  ${BOOST_INCLUDE_DIRS}
  ${GEANT4_INCLUDE_DIR}
)

add_subdirectory(geometries)
add_subdirectory(materials)
add_subdirectory(processes)
add_subdirectory(physics_lists)
add_subdirectory(primaries)
add_subdirectory(utils)
add_subdirectory(python)
