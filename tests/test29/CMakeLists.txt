cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test29)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test29 test29.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test29 COMMAND test29 ${CMAKE_CURRENT_SOURCE_DIR}/test29.in
                       BUILD test29
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
