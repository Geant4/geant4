cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test65)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})
##include(UseGeant4_internal)

GEANT4_EXECUTABLE(test65 test65.cc src/*.cc)
GEANT4_EXECUTABLE(test65_hadronic_exerciser test65_hadronic_exerciser.cc)

add_dependencies(test65 test65_hadronic_exerciser) # Artificial to force serialization (problem on windows)

if(NOT GEANT4_TEST_ENVIRONMENT MATCHES G4LENDDATA) 
    message(STATUS "G4 TESTS: G4LENDDATA ennvironment var not set. --> test65 disabled")  
  return()
endif()
#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test65_hadronic_exerciser COMMAND test65_hadronic_exerciser
                                          OUTPUT  test65.in
                                          BUILD test65_hadronic_exerciser)
GEANT4_ADD_TEST(test65 COMMAND test65 test65.in
                       BUILD test65
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
                       DEPENDS test65_hadronic_exerciser)
