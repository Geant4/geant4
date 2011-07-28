cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test37)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test37 test37.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test37 COMMAND test37 ${CMAKE_CURRENT_SOURCE_DIR}/test37.in
                       BUILD test37
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
