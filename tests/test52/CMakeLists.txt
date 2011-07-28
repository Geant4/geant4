cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test52)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

find_package(AIDA QUIET)

if(NOT AIDA_FOUND)
  message(STATUS "G4 TESTS: AIDA package not found. --> test52 disabled")  
  return()
endif()

include_directories(${AIDA_INCLUDE_DIR})
GEANT4_EXECUTABLE(test52 test52.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test52 COMMAND test52 ${CMAKE_CURRENT_SOURCE_DIR}/standard.mac
                       BUILD test52
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
