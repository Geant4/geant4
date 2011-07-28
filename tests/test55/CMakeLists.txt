cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test55)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test55 test55.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test55 COMMAND test55 ${CMAKE_CURRENT_SOURCE_DIR}/test55.in
                       BUILD test55
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
