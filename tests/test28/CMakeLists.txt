cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test28)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test28 test28.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test28 COMMAND test28 ${CMAKE_CURRENT_SOURCE_DIR}/test28.in
                       BUILD test28
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
