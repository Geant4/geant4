# - CmakeLists.txt for parallel examples

cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

find_package(Geant4)
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------

find_package(MPI QUIET)
if(MPI_CXX_FOUND)
  add_subdirectory(MPI)
else()  
  message(STATUS "G4 Examples: MPI not installed --> MPI examples disabled")  
endif()

find_package(TBB QUIET)
if(TBB_FOUND)
  add_subdirectory(TBB)
else()  
  message(STATUS "G4 Examples: TBB not installed --> TBB example disabled")  
endif()

