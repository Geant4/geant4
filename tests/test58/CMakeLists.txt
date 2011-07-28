cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test58)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test58 test58.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test58 COMMAND test58 ${CMAKE_CURRENT_SOURCE_DIR}/test58.in
                       BUILD test58
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
