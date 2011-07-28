cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test22)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test22 test22.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test22 COMMAND test22 ${CMAKE_CURRENT_SOURCE_DIR}/test22.in
                       BUILD test22
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
GEANT4_ADD_TEST(test22-largeN COMMAND test22 ${CMAKE_CURRENT_SOURCE_DIR}/test22.large_N.in
                       DEPENDS test22
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
