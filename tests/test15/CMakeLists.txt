cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test15)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test15 test15.cc src/*.cc)
GEANT4_EXECUTABLE(test15_hadronic_exerciser test15_hadronic_exerciser.cc)
add_dependencies(test15 test15_hadronic_exerciser) # Artificial to force serialization (problem on windows)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test15-build BUILD test15)
GEANT4_ADD_TEST(test15_hadronic_exerciser COMMAND test15_hadronic_exerciser
                                          OUTPUT  test15.in
                                          DEPENDS test15-build)
GEANT4_ADD_TEST(test15 COMMAND test15 test15.in
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
                       DEPENDS test15_hadronic_exerciser)
