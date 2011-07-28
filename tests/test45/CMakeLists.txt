cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test45)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test45 test45.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test45 COMMAND test45 ${CMAKE_CURRENT_SOURCE_DIR}/test45.in
                       BUILD test45
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
