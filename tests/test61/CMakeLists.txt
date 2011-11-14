cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test61)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test61 test61.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test61 COMMAND test61 ${CMAKE_CURRENT_SOURCE_DIR}/test61.in
                       BUILD test61 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

