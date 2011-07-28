cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test44)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test44 test44.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test44 COMMAND test44 ${CMAKE_CURRENT_SOURCE_DIR}/test44.in
                       BUILD test44
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
