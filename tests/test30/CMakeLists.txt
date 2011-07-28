cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test30)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test30 test30.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test30 COMMAND test30 ${CMAKE_CURRENT_SOURCE_DIR}/test30.in
                       BUILD test30
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
