cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test05)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test05 test05.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test05 COMMAND test05 ${CMAKE_CURRENT_SOURCE_DIR}/test05.in
                       BUILD test05
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
GEANT4_ADD_TEST(test05-largeN COMMAND test05 ${CMAKE_CURRENT_SOURCE_DIR}/test05.largeN.in
                       DEPENDS test05
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
