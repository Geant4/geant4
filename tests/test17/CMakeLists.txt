cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test17)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test17 test17.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test17 COMMAND test17 ${CMAKE_CURRENT_SOURCE_DIR}/test17.in
                       BUILD test17
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
GEANT4_ADD_TEST(test17-largeN COMMAND test17 ${CMAKE_CURRENT_SOURCE_DIR}/test17.large_N.in
                       DEPENDS test17
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
