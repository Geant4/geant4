# - CmakeLists.txt for MPI examples
cmake_minimum_required(VERSION 3.16...3.27)
project(G4mpi)

# - Require both MPI and Geant4.
# Only MPI component required is CXX (ability to use C++ with C API).
find_package(MPI REQUIRED CXX)
find_package(Geant4 REQUIRED ui_all vis_all)

# - The G4mpi library used by all the examples
add_subdirectory(source)

# - The example applications themselves
add_subdirectory(exMPI01)
add_subdirectory(exMPI02)
add_subdirectory(exMPI03)
add_subdirectory(exMPI04)


