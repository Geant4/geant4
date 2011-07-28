cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test36)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test36 test36.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test36 COMMAND test36 ${CMAKE_CURRENT_SOURCE_DIR}/test36.in
                       BUILD test36
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
