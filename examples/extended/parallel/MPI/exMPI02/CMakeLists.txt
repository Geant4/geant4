#----------------------------------------------------------------------------
# Create a new project level to avoid name clashes with other exMPI0N
project(exMPI02)

# Find ROOT and return if we don't have it to prvent building this project
find_package(ROOT QUIET COMPONENTS Hist RIO)
if(NOT ROOT_FOUND)
  message(STATUS "ROOT not found, disabling exMPI02 build")
  return()
endif()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(exMPI02 exMPI02.cc ${sources} ${headers})
target_include_directories(exMPI02 PRIVATE include)
target_link_libraries(exMPI02 PRIVATE ${Geant4_LIBRARIES} G4mpi::G4mpi ROOT::Hist ROOT::RIO)

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build exMPI02. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(exMPI02_SCRIPTS
    run.mac
    vis.mac
  )

foreach(_script ${exMPI02_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

