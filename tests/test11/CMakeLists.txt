cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test11)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test11 test11.cc src/*.cc)
GEANT4_EXECUTABLE(test11_hadronic_exerciser test11_hadronic_exerciser.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test11_hadronic_exerciser COMMAND test11_hadronic_exerciser
                                          OUTPUT  test11.in
                                          BUILD test11_hadronic_exerciser)
GEANT4_ADD_TEST(test11 COMMAND test11 test11.in
                       BUILD test11
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
                       DEPENDS test11_hadronic_exerciser)



