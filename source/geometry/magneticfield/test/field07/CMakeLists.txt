cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(field07)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})
##include(UseGeant4_internal)

GEANT4_EXECUTABLE(field07 field07.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(field07-build BUILD test11)
GEANT4_ADD_TEST(field07 COMMAND field07 field07.in
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT} )



