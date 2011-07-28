cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test13)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test13 test13.cc src/*.cc)
GEANT4_EXECUTABLE(test13_hadronic_exerciser test13_hadronic_exerciser.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test13_hadronic_exerciser COMMAND test13_hadronic_exerciser
                                          OUTPUT  test13.in
                                          BUILD test13_hadronic_exerciser)
GEANT4_ADD_TEST(test13 COMMAND test13 test13.in
                       BUILD test13
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
                       DEPENDS test13_hadronic_exerciser)
