cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test34)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test34 test34.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test34 COMMAND test34 ${CMAKE_CURRENT_SOURCE_DIR}/test34.in
                       BUILD test34
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
