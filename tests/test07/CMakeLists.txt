cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test07)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test07 test07.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test07 COMMAND test07 ${CMAKE_CURRENT_SOURCE_DIR}/test07.in
                       BUILD test07
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
GEANT4_ADD_TEST(test07-largeN COMMAND test07 ${CMAKE_CURRENT_SOURCE_DIR}/test07.large_N.in
                       DEPENDS test07
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

