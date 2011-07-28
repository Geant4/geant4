cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test40)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test40 test40.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test40 COMMAND test40 ${CMAKE_CURRENT_SOURCE_DIR}/test40.in
                       BUILD test40
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
GEANT4_ADD_TEST(test40-largeN COMMAND test40 ${CMAKE_CURRENT_SOURCE_DIR}/test40.large_N.in
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
