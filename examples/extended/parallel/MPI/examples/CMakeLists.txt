# - CmakeLists.txt for MPI examples

cmake_minimum_required(VERSION 3.8...3.18)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

add_subdirectory(exMPI01)
add_subdirectory(exMPI02)
add_subdirectory(exMPI03)
add_subdirectory(exMPI04)
