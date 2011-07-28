cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test46)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test46 test46.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test46 COMMAND test46 ${CMAKE_CURRENT_SOURCE_DIR}/test46.in
                       BUILD test46
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
