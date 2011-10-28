cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test69)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test69 test69.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test69 COMMAND test69 ${CMAKE_CURRENT_SOURCE_DIR}/test69.in
                       BUILD test69
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

