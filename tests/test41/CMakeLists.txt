cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test41)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test41 test41.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test41 COMMAND test41 ${CMAKE_CURRENT_SOURCE_DIR}/test41.in
                       BUILD test41
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
