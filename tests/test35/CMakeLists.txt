cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test35)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test35 test35.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test35 COMMAND test35 ${CMAKE_CURRENT_SOURCE_DIR}/test35.in
                       BUILD test35
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
