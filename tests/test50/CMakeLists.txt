cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test50)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test50 test50.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test50 COMMAND test50 ${CMAKE_CURRENT_SOURCE_DIR}/test50.in
                       BUILD test50
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
