cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test15)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test15 test15.cc src/*.cc)
GEANT4_EXECUTABLE(test15_hadronic_exerciser test15_hadronic_exerciser.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test15_hadronic_exerciser COMMAND test15_hadronic_exerciser
                                          OUTPUT  test15.in
                                          BUILD test15_hadronic_exerciser)
GEANT4_ADD_TEST(test15 COMMAND test15 test15.in
                       BUILD test15
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
                       DEPENDS test15_hadronic_exerciser)
