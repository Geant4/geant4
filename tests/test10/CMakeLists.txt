cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test10)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test10 test10.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test10 COMMAND test10 ${CMAKE_CURRENT_SOURCE_DIR}/test10.in
                       BUILD test10
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
GEANT4_ADD_TEST(test10-largeN COMMAND test10 ${CMAKE_CURRENT_SOURCE_DIR}/test10.large_N.in
                       DEPENDS test10
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

