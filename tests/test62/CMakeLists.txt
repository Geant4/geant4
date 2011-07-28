cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test62)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test62 test62.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test62 COMMAND test62 ${CMAKE_CURRENT_SOURCE_DIR}/test62.in
                       BUILD test62
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

