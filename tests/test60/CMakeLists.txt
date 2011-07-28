cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test60)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test60 test60.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test60 COMMAND test60 ${CMAKE_CURRENT_SOURCE_DIR}/test60.in
                       BUILD test60
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
GEANT4_ADD_TEST(test60-largeN COMMAND test60 ${CMAKE_CURRENT_SOURCE_DIR}/test60.large_N.in
                       BUILD test60
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

