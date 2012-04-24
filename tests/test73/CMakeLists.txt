cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test73)

find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

find_package(ROOT QUIET)
if(NOT ROOT_FOUND)
  message(STATUS "G4 TESTS: ROOT package not found. --> test73 disabled")
  return()
endif()


add_definitions(-DG4ANALYSIS_USEROOT)
include_directories(${ROOT_INCLUDE_DIR})
GEANT4_EXECUTABLE(test73 test73.cc src/*.cc LIBRARIES ${ROOT_LIBRARIES})

#---Test definition-----------------------------------------
file(GLOB MACROS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.mac xmacro.C)
file(COPY ${MACROS} DESTINATION .)

GEANT4_ADD_TEST(test73-FTFP_BERT COMMAND test73 ${CMAKE_CURRENT_SOURCE_DIR}/run_all.mac FTFP_BERT
		       BUILD test73
		       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
