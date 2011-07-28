cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test19)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test19 test19.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test19 COMMAND test19 ${CMAKE_CURRENT_SOURCE_DIR}/chipstest.in
                       BUILD test19
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
