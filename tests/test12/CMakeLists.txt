cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test12)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test12 test12.cc src/*.cc)
GEANT4_EXECUTABLE(test12_hadronic_exerciser test12_hadronic_exerciser.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test12_hadronic_exerciser COMMAND test12_hadronic_exerciser
                                          OUTPUT  test12.in
                                          BUILD test12_hadronic_exerciser)
GEANT4_ADD_TEST(test12 COMMAND test12 test12.in
                       BUILD test12
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
                       DEPENDS test12_hadronic_exerciser)
