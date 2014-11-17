# - add libs components

set(G4MODULES_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR}/Geant4)

include_directories (
  ${PYTHON_INCLUDE_PATH}
  ${Boost_INCLUDE_DIRS}
  ${PROJECT_SOURCE_DIR}/source/boost
  ${GEANT4_INCLUDE_DIR}
)

link_directories (${GEANT4_LIBRARY_DIR} ${Boost_LIBRARY_DIRS})

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
add_subdirectory(gdml)

if (${PYTHON_VERSION_MAJOR} MATCHES "2")
  add_subdirectory(python)
elseif (${PYTHON_VERSION_MAJOR} MATCHES "3")
  add_subdirectory(python3)
endif()
