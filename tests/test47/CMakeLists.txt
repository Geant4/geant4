cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test47)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

find_package(ROOT QUIET)
if(NOT ROOT_FOUND)
  message(STATUS "G4 TESTS: ROOT package not found. --> test47 disabled")  
  return()
endif()

include_directories(${ROOT_INCLUDE_DIR})
GEANT4_EXECUTABLE(test47 test47.cc src/*.cc LIBRARIES ${ROOT_LIBRARIES})

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test47 COMMAND test47 ${CMAKE_CURRENT_SOURCE_DIR}/test47.inp
                       BUILD test47
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT} ROOTSYS=${ROOTSYS})
