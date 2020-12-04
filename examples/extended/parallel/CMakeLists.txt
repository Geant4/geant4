# - CmakeLists.txt for parallel examples

cmake_minimum_required(VERSION 3.8...3.18)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

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

