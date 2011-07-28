cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test21)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test21 test21.cc src/*.cc)
GEANT4_EXECUTABLE(test21_hadronic_exerciser test21_hadronic_exerciser.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test21_hadronic_exerciser COMMAND test21_hadronic_exerciser
                                          OUTPUT  test21.in
                                          BUILD test21_hadronic_exerciser)
GEANT4_ADD_TEST(test21 COMMAND test21 test21.in
                       BUILD test21
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
                       DEPENDS test21_hadronic_exerciser)
