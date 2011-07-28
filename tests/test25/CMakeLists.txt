cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test25)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test25 test25.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test25 COMMAND test25 ${CMAKE_CURRENT_SOURCE_DIR}/test25.in
                       BUILD test25
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
