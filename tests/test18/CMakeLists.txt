cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test18)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test18 test18.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test18 COMMAND test18 ${CMAKE_CURRENT_SOURCE_DIR}/test18.in
                       BUILD test18
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
GEANT4_ADD_TEST(test18-largeN COMMAND test18 ${CMAKE_CURRENT_SOURCE_DIR}/test18.large_N.in
                       DEPENDS test18
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
