#----------------------------------------------------------------------------
# Create a new project level to avoid name clashes with other exMPI0N
project(exMPI04)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(exMPI04 exMPI04.cc ${sources} ${headers})
target_include_directories(exMPI04 PRIVATE include)
target_link_libraries(exMPI04 PRIVATE ${Geant4_LIBRARIES} G4mpi::G4mpi)

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build exMPI04. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(exMPI04_SCRIPTS
    run.mac
    vis.mac
  )

foreach(_script ${exMPI04_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

