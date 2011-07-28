cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test67)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test67 test67.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test67 COMMAND test67 ${CMAKE_CURRENT_SOURCE_DIR}/test67.in
                       BUILD test67
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
