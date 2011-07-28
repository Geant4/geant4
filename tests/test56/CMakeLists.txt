cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test56)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test56 test56.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test56 COMMAND test56 ${CMAKE_CURRENT_SOURCE_DIR}/test56.in
                       BUILD test56
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
