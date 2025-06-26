# - CmakeLists.txt for parallel examples

cmake_minimum_required(VERSION 3.16...3.27)

find_package(Geant4)
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------

find_package(MPI QUIET)
if(MPI_CXX_FOUND)
  add_subdirectory(MPI)
else()  
  message(STATUS "G4 Examples: MPI not installed --> MPI examples disabled")  
endif()


