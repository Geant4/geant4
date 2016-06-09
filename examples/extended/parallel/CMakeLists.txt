# - CmakeLists.txt for parallel examples

cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

find_package(MPI QUIET)
if(MPI_CXX_FOUND)
  add_subdirectory(MPI)
else()  
  message(STATUS "G4 Examples: MPI not installed --> MPI examples disabled")  
endif()

