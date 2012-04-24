cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test74)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})
GEANT4_EXECUTABLE(test74 test74.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test74 COMMAND test74 ${CMAKE_CURRENT_SOURCE_DIR}/test74.in
                       BUILD test74
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
GEANT4_ADD_TEST(test74-largeN COMMAND test74 ${CMAKE_CURRENT_SOURCE_DIR}/test74.large_N.in
                       DEPENDS test74
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

