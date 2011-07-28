cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test09)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test09 test09.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test09 COMMAND test09 ${CMAKE_CURRENT_SOURCE_DIR}/test09.in
                       BUILD test09
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
GEANT4_ADD_TEST(test09-largeN COMMAND test09 ${CMAKE_CURRENT_SOURCE_DIR}/test09.large_N.in
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
