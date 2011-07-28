cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test31)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test31 test31.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test31 COMMAND test31 ${CMAKE_CURRENT_SOURCE_DIR}/test31.in
                       BUILD test31
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
