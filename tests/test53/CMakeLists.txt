cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test53)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test53 test53.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test53 COMMAND test53 ${CMAKE_CURRENT_SOURCE_DIR}/test53.in
                       BUILD test53
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
