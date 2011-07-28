cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test27)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test27 test27.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test27 COMMAND test27 ${CMAKE_CURRENT_SOURCE_DIR}/test27.in
                       BUILD test27
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
