cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test39)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test39 test39.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test39 COMMAND test39 ${CMAKE_CURRENT_SOURCE_DIR}/test39.in
                       BUILD test39
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
