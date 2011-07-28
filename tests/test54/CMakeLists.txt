cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test54)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test54 test54.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test54 COMMAND test54 ${CMAKE_CURRENT_SOURCE_DIR}/test54.in
                       BUILD test54
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
