cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(exampleB02)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

find_package(AIDA QUIET)

if(NOT AIDA_FOUND)
  message(STATUS "G4 TESTS: AIDA package not found. --> exampleB02 disabled")  
  return()
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${Geant4_INCLUDE_DIR})
file(GLOB sources ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc)

add_executable(exampleB02 EXCLUDE_FROM_ALL exampleB02.cc ${sources})
target_link_libraries(exampleB02 ${Geant4_LIBRARIES})

