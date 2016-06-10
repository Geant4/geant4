# - CmakeLists.txt for MPI examples
# This file should be used when configuring Geant4 with GEANT4_BUILD_EXAMPLES=ON
# Do not use in other cases
# To compile G4mpi as a standalone package, first compile the source subdir
# and then examples, see included README files
project(G4mpi)

cmake_minimum_required(VERSION 2.8.3 FATAL_ERROR)

add_subdirectory(source)
set(G4mpi_DIR "${PROJECT_BINARY_DIR}/source")
add_subdirectory(examples)
add_dependencies(exMPI01 G4mpi)
if(TARGET exMPI02) #If ROOT not found, this does not exists
   add_dependencies(exMPI02 G4mpi)
endif()
add_dependencies(exMPI03 G4mpi)


add_custom_target(MPI DEPENDS G4mpi exMPI01 exMPI02 exMPI03)


