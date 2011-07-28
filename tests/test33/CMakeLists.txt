cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test33)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test33 test33.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
configure_file(SamplingChange.mac SamplingChange.mac COPYONLY) # test expects file to be local
GEANT4_ADD_TEST(test33 COMMAND test33 ${CMAKE_CURRENT_SOURCE_DIR}/MassGeo_TimedApp.mac
                       BUILD test33
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
GEANT4_ADD_TEST(test33-1 COMMAND test33 ${CMAKE_CURRENT_SOURCE_DIR}/ParallelGeo_TimedApp.mac
                       DEPENDS test33
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
