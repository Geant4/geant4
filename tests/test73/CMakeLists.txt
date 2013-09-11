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
#Prepare running are: 
# 1- create working directories
# 2- Copy macros in the directories
# 3- Define tests (since we run in subdirs we can run in parallel)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/local ${CMAKE_CURRENT_BINARY_DIR}/FTFP_BERT_EMX ${CMAKE_CURRENT_BINARY_DIR}/FTFP_BERT)
file(GLOB MACROS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.mac xmacro.C)
file(COPY ${MACROS} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/local)
file(COPY ${MACROS} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/FTFP_BERT_EMX)
file(COPY ${MACROS} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/FTFP_BERT)

GEANT4_ADD_TEST( test73-build BUILD test73 LABELS Nightly Continuous PhysicsChecks) 

#The following tests are enabled only in PhysicsChecks since they 
#Will fail for physics reasons: the MSC model is not adeguate for pi-
#When new G4 developments are such that they will not fail anymore
#they can be promoted to nightlies adding the Nichtly and Continous
#LABELS

#GEANT4_ADD_TEST(test73-local COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test73 run_all.mac local
#	       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/local
#	       DEPENDS test73-build
#	       LABELS PhysicsChecks
#	       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})


GEANT4_ADD_TEST(test73-FTFP_BERT_EMX COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test73 run_all.mac FTFP_BERT_EMX
	       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/FTFP_BERT_EMX
	       DEPENDS test73-build
	       LABELS PhysicsChecks
	       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

GEANT4_ADD_TEST(test73-FTFP_BERT COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test73 run_all.mac FTFP_BERT
	       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/FTFP_BERT
	       DEPENDS test73-build
	       LABELS PhysicsChecks
	       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

