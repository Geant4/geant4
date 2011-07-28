cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test49)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test49 test49.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test49 COMMAND test49 ${CMAKE_CURRENT_SOURCE_DIR}/chipstest.in
                       BUILD test49
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
