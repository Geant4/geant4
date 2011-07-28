cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test24)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test24 test24.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test24 COMMAND test24 ${CMAKE_CURRENT_SOURCE_DIR}/test24.in
                       BUILD test24
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
